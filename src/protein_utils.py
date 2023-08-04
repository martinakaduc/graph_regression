import numpy as np
import requests
import os

from tqdm import tqdm
from collections import defaultdict
from rdkit.Chem import DetectChemistryProblems
from rdkit.Chem.rdmolfiles import MolFromPDBFile, MolFromSmiles

import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
from requests.adapters import HTTPAdapter, Retry

from constants import PDB_FOLDER, DATA_FOLDER, UNIPROT_API_URL, POLLING_INTERVAL
from utils import readCSV, findMatchingData

session = requests.Session()
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session.mount("https://", HTTPAdapter(max_retries=retries))

def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise
        
def get_id_mapping_results_link(job_id):
    url = f"{UNIPROT_API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{UNIPROT_API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]

def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results

def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)
        
def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    # print(f"Fetched: {n_fetched} / {total}")
    
def get_id_mapping_results_search(url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results

def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{UNIPROT_API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                # print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])
        
# Function for data processing
def getProteinFASTA(uniprot_id):
    baseUrl = "http://www.uniprot.org/uniprot/"
    currentUrl = baseUrl + uniprot_id + ".fasta"
    response = requests.post(currentUrl)
    fasta = response.text
    return fasta

def getProteinRCSB(pdb_id, uniprot_id, overwrite=False):
    url = f"http://files.rcsb.org/view/{pdb_id}.pdb"
    filename = f'{uniprot_id}.pdb'
    cache_path = os.path.join(PDB_FOLDER, filename)
    if os.path.exists(cache_path) and not overwrite:
        return cache_path

    pdb_req = requests.get(url)
    try:
        pdb_req.raise_for_status()
    except requests.exceptions.HTTPError:
        return None
    
    open(cache_path, 'w').write(pdb_req.text)
    return cache_path

def converUniProt2PDB(uniprot_ids):
    '''
    curl --form 'from="UniProtKB_AC-ID"' \
     --form 'to="PDB"' \
     --form 'ids="P07550"' \
     https://rest.uniprot.org/idmapping/run
    '''
    if isinstance(uniprot_ids, str):
        uniprot_ids = [uniprot_ids]
        
    job_id = submit_id_mapping(
        from_db="UniProtKB_AC-ID", to_db="PDB", ids=uniprot_ids
    )
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
        
    return results['results']
    
def getProteinMol(uniprot_id, pdb_ids, overwrite=False):
#     protein_fasta = getProteinFASTA(uniprot_id)
#     protein_mol = MolFromFASTA(protein_fasta)
#     MolToPDBFile(protein_mol, uniprot_id + ".pdb")

    if not overwrite:
        protein_path = getProteinRCSB(0, uniprot_id, overwrite=overwrite)
        protein_mol = MolFromPDBFile(protein_path)
        return protein_mol
        
    # Trying to find valid protein structure
    # pdb_ids = converUniProt2PDB(uniprot_id)
    # assert pdb_ids is not None
    if len(pdb_ids) == 0:
        return None
    
    id_index = 0
    pdb_id = pdb_ids[id_index]
    protein_path = getProteinRCSB(pdb_id, uniprot_id, overwrite=overwrite)
    if protein_path is None:
        protein_mol = None
    else:
        protein_mol = MolFromPDBFile(protein_path)
        if protein_mol is not None and len(DetectChemistryProblems(protein_mol)) > 0:
            # Ensure no problem in protein PDB
            protein_mol = None
    
    # Loop until found a valid structure
    while protein_mol is None and id_index + 1 < len(pdb_ids):
        id_index += 1
        pdb_id = pdb_ids[id_index]
        protein_path = getProteinRCSB(pdb_id, uniprot_id, overwrite=overwrite)
        if protein_path is None:
            protein_mol = None
        else:
            protein_mol = MolFromPDBFile(protein_path)
        
        # Ensure no problem in protein PDB
        if protein_mol is not None and len(DetectChemistryProblems(protein_mol)) > 0:
            protein_mol = None
            
    if protein_mol is None:
        filename = f'{uniprot_id}.pdb'
        cache_path = os.path.join(PDB_FOLDER, filename)
        os.system("rm -rf %s" % cache_path)
        
    return protein_mol

def getLigandMol(ligand_smile):
    return MolFromSmiles(ligand_smile)

def convertProtein2UniprotID(mapping_df, protein_ids):
    converted_uniprot_ids = []
    for protein_id in protein_ids:
        uniprot_ids = findMatchingData(mapping_df, protein_id, "GtoPdb IUPHAR ID", "UniProtKB ID").tolist()
        if len(uniprot_ids) > 0:
            # Just pick the first protein variant
            converted_uniprot_ids.append(uniprot_ids[0])
    return converted_uniprot_ids

# Get protein by family and spliting training/testing
def getProteinsByFamilies(family_target_mapping_file):
    family_target_df = readCSV(family_target_mapping_file)
    families = set(family_target_df["Family id"].tolist())
    
    GtP_to_UniProt_mapping_file = os.path.join(DATA_FOLDER, "GtP_to_UniProt_mapping.csv")
    GtP_to_UniProt_mapping_df = readCSV(GtP_to_UniProt_mapping_file)
    
    proteins_by_families = defaultdict(list)
    protein_family_mapping = defaultdict(str)
    
    print("Extracting protein families...")
    for family in tqdm(families):
        proteins_ids = findMatchingData(family_target_df, family,"Family id", "Target id").tolist()
        proteins_by_families[family] = convertProtein2UniprotID(GtP_to_UniProt_mapping_df, proteins_ids)
        
        for pid in proteins_by_families[family]:
            protein_family_mapping[pid] = family
            
    # Convert UniProt ID to PDB ID
    print("Converting UniProt ID to PDB ID...")
    all_protein_uniprot_ids = list(protein_family_mapping.keys())
    uniprot2pdb_mapping = converUniProt2PDB(all_protein_uniprot_ids)
    uniprot2pdb_mapping_dict = defaultdict(list)
    for mapping in uniprot2pdb_mapping:
        uniprot_id, pdb_id = mapping['from'], mapping['to']
        uniprot2pdb_mapping_dict[uniprot_id].append(pdb_id)
    
    return proteins_by_families, protein_family_mapping, uniprot2pdb_mapping_dict


def getTrainTestByProteins(proteins_by_families, uniprot2pdb_mapping_dict, test_percentage=0.33):
    training_set = set([])
    testing_set = set([])
        
    print("Fetching proteins by families from PDB...")
    for family, protein_uniprot_ids in tqdm(proteins_by_families.items()):
        if len(protein_uniprot_ids) < 2:
            continue
            
        print("Protein family:", family)
        print("Number of proteins:", len(protein_uniprot_ids))
        
        p_test_idx = 0
        test_0 = None
        while test_0 is None and p_test_idx < len(protein_uniprot_ids):
            pdb_ids = uniprot2pdb_mapping_dict[protein_uniprot_ids[p_test_idx]]
            test_0 = getProteinMol(protein_uniprot_ids[p_test_idx], pdb_ids, overwrite=True)
            p_test_idx += 1

        if p_test_idx >= len(protein_uniprot_ids):
            # No protein left for training...
            # So we ignore this family
            continue

        p_train_idx = p_test_idx
        train_0 = None
        while train_0 is None and p_train_idx < len(protein_uniprot_ids):
            pdb_ids = uniprot2pdb_mapping_dict[protein_uniprot_ids[p_train_idx]]
            train_0 = getProteinMol(protein_uniprot_ids[p_train_idx], pdb_ids, overwrite=True)
            p_train_idx += 1

        if train_0 is None:
            # No protein valid for training
            continue

        # If we reach here, there are at least 1 valid protein for training and 1 valid one for testing
        p_test_idx = p_test_idx - 1
        p_train_idx = p_train_idx -1
        training_set.add(protein_uniprot_ids[p_train_idx])
        testing_set.add(protein_uniprot_ids[p_test_idx])

        # The rest proteins will be randomly put into training/testing
        p_idx = p_train_idx + 1
        while p_idx < len(protein_uniprot_ids):
            pdb_ids = uniprot2pdb_mapping_dict[protein_uniprot_ids[p_idx]]
            valid_protein = getProteinMol(protein_uniprot_ids[p_idx], pdb_ids, overwrite=True)
            
            if valid_protein:
                into_test_set = np.random.choice([1,0], p=[test_percentage, 1-test_percentage])
                if into_test_set:
                    testing_set.add(protein_uniprot_ids[p_idx])
                else:
                    training_set.add(protein_uniprot_ids[p_idx])
            p_idx += 1
            
    return list(training_set), list(testing_set)