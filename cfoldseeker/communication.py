#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import requests
import time
import sys
from pathlib import Path
from tqdm.contrib.concurrent import thread_map

LOG = logging.getLogger(__name__)


"""
Submits one structure file to the FoldSeek API and returns the submission ticket in dictionary form.
"""
def submit_foldseek_query(query_path: Path, dbs: list, taxfilters: list) -> None:
    FOLDSEEK_SUBMISSION_URL = "https://search.foldseek.com/api/ticket"
    
    with open(query_path, "rb") as f:
        files = {"q": f}
        data = [("mode", "3diaa")]
        for db in dbs:
            data.append(("database[]", db))
        for taxfilt in taxfilters:
            data.append(('taxfilter', taxfilt))
        LOG.debug(f"Posting request on FoldSeek webserver {FOLDSEEK_SUBMISSION_URL}")
        LOG.debug(f"using the following parameters: {data}")
        response = requests.post(FOLDSEEK_SUBMISSION_URL, files=files, data=data)
        if response.status_code == 200:
            LOG.debug(f'Query {query_path} successfully submitted!')
            return response.json()
        else:
            LOG.critical(f"Error submitting query {query_path}!")
            sys.exit()
            
    return None


"""
Checks the status of the job associated with the provided job ID. Returns the status flag.
"""
def check_query_status(job_id: str) -> str:
    FOLDSEEK_SUBMISSION_URL = "https://search.foldseek.com/api/ticket"
    
    url = f"{FOLDSEEK_SUBMISSION_URL}/{job_id}"
    LOG.debug(f'Checking out status of job {job_id}')
    results = requests.get(url).json()
    status = results['status']
    
    return status


"""
Waits for the FoldSeek job associated with the provided job ID to complete, and retrieves its result.
Returns the parsed results in a dictionary.
"""
def retrieve_foldseek_results(job_id: str) -> dict:
    FOLDSEEK_RESULTS_URL = "https://search.foldseek.com/api/result"
    
    while True:
        status = check_query_status(job_id)
        if status == "COMPLETE":
            LOG.debug(f"Job {job_id} has completed! Downloading results...")
            entry = 0
            url = f"{FOLDSEEK_RESULTS_URL}/{job_id}/{entry}"
            results = requests.get(url).json()
            break
        else:
            LOG.debug(f"Job {job_id} has not completed yet. Waiting another 10 seconds...")
            time.sleep(10)
            
    return results


def pull_from_ena(entry: str, max_retries: int = 3) -> None | str:
    """
    Pulls a GenPept record from the ENA Browser API.
    """
    ENA_BROWSER_URL = "https://www.ebi.ac.uk/ena/browser/api/embl"
    trials = 0
    
    url = f"{ENA_BROWSER_URL}/{entry}"
    LOG.debug(f'Going to pull GenPept record from {url}')
    while trials < max_retries:
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        elif response.status_code == 429:
            trials += 1
            time.sleep(2)
        else:
            LOG.warning(f'Error pulling GenPept entry {entry}, code returned: {response.status_code}')
            return None


def pull_from_unisave(entry: str, max_retries: int = 3) -> None | str:
    """
    Pulls a UniSave record from UniProt REST API.
    """
    UNISAVE_REST_URL = "https://rest.uniprot.org/unisave"
    trials = 0
    
    url = f"{UNISAVE_REST_URL}/{entry}?format=txt"
    LOG.debug(f'Going to pull UniSave entry from {url}')
    while trials < max_retries:
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        elif response.status_code == 429:
            trials += 1
            time.sleep(2)
        else:
            LOG.warning(f'Error pulling UniSave record {entry}, code returned: {response.status_code}')
            return None


def pull_dict_from_unisave(entries: list, max_workers: int = 1, no_progress: bool = False) -> dict:
    """
    Pulls a list of UniSave entries and returns it in a dictionary.
    """
    LOG.info(f'Going to pull {len(entries)} UniSave records')
    
    unisave_entries = thread_map(pull_from_unisave, entries, 
                                 max_workers = max_workers,
                                 leave = False,
                                 disable = no_progress)
    unisave_entries = [e for e in unisave_entries if e != None]
    return dict(zip(entries, unisave_entries))

