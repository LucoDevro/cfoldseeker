#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import requests
import time
from pathlib import Path
from tqdm.contrib.concurrent import thread_map
from tqdm.contrib.logging import logging_redirect_tqdm


LOG = logging.getLogger(__name__)


def submit_foldseek_query(query_path: Path, dbs: list, taxfilters: list) -> dict:
    """
    Submits a structure file to the FoldSeek API for processing.
    
    Sends a protein structure query to the FoldSeek webserver with specified
    databases and taxonomic filters. Returns the submission ticket on success
    or raises an error on failure.
    
    Args:
        query_path: Path object pointing to the structure file to submit.
        dbs: List of database names to search against.
        taxfilters: List of taxonomic filters to apply to the search.
    
    Returns:
        A dictionary containing the submission ticket and metadata from the
        FoldSeek API response.
        
    Raises:
        RuntimeError: If submission fails.
        
    """
    FOLDSEEK_SUBMISSION_URL = "https://search.foldseek.com/api/ticket"
    
    with open(query_path, "rb") as f:
        # Compose ticket
        files = {"q": f}
        data = [("mode", "3diaa")]
        for db in dbs:
            data.append(("database[]", db))
        for taxfilt in taxfilters:
            data.append(('taxfilter', taxfilt))
        
        # Submit ticket
        LOG.debug(f"Posting request on FoldSeek webserver {FOLDSEEK_SUBMISSION_URL}")
        LOG.debug(f"using the following parameters: {data}")
        response = requests.post(FOLDSEEK_SUBMISSION_URL, files=files, data=data)
        
        # Checking server response
        if response.status_code == 200:
            LOG.debug(f'Query {query_path} successfully submitted!')
            return response.json()
        else:
            msg = f"Error submitting query {query_path}! Status code: {response.status_code}"
            LOG.critical(msg)
            raise RuntimeError(msg)
            
    return None


def check_query_status(job_id: str) -> str:
    """
    Retrieves the current status of a FoldSeek job.
    
    Queries the FoldSeek API to check the processing status of a previously
    submitted job using its unique job ID.
    
    Args:
        job_id: The unique identifier for the FoldSeek job.
    
    Returns:
        A string indicating the job status (e.g., "COMPLETE", "RUNNING", etc.).
    """
    FOLDSEEK_SUBMISSION_URL = "https://search.foldseek.com/api/ticket"
    
    url = f"{FOLDSEEK_SUBMISSION_URL}/{job_id}"
    LOG.debug(f'Checking out status of job {job_id}')
    results = requests.get(url).json()
    status = results['status']
    
    return status


def retrieve_foldseek_results(job_id: str) -> dict:
    """
    Waits for a FoldSeek job to complete and retrieves its results.
    
    Polls the job status at regular intervals until completion, then downloads
    and returns the parsed results from the FoldSeek API.
    
    Args:
        job_id: The unique identifier for the FoldSeek job.
    
    Returns:
        A dictionary containing the parsed results from the completed FoldSeek job.
    """
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
    Retrieves a GenPept record from the ENA Browser API.
    
    Attempts to fetch a GenPept sequence record from the European Nucleotide
    Archive (ENA) with retry logic for rate-limited responses.
    
    Args:
        entry: The accession number or identifier of the GenPept record to retrieve.
        max_retries: Maximum number of retry attempts for rate-limited requests (error code 429).
            Defaults to 3. Waiting time between trials is 5 seconds.
    
    Returns:
        A string containing the GenPept record in text format, or None if the
        retrieval fails after max retries or an unexpected error occurs.
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
            time.sleep(5)
        else:
            LOG.warning(f'Error pulling GenPept entry {entry}, code returned: {response.status_code}')
            return None


def pull_from_unisave(entry: str, max_retries: int = 3) -> None | str:
    """
    Retrieves a UniSave record from the UniProt REST API.
    
    Fetches a protein sequence record from UniSave (UniProt archive) with
    retry logic for rate-limited responses.
    
    Args:
        entry: The UniProt accession number of the record to retrieve.
        max_retries: Maximum number of retry attempts for rate-limited requests (429).
            Defaults to 3. Waiting time between trials is 5 seconds.
    
    Returns:
        A string containing the UniSave record in text format, or None if the
        retrieval fails after max retries or an unexpected error occurs.
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
            time.sleep(5)
        else:
            LOG.warning(f'Error pulling UniSave record {entry}, code returned: {response.status_code}')
            return None


def pull_dict_from_unisave(entries: list, max_workers: int = 1, no_progress: bool = False) -> dict:
    """
    Retrieves multiple UniSave records and returns them as a dictionary.
    
    Fetches a list of UniSave entries in parallel and returns them mapped to
    their original accession numbers. Failed retrievals are filtered out.
    
    Args:
        entries: List of UniProt accession numbers to retrieve.
        max_workers: Number of worker threads for parallel retrieval.
            Defaults to 1.
        no_progress: If True, suppresses the progress bar during retrieval.
            Defaults to False.
    
    Returns:
        A dictionary mapping each successfully retrieved accession number to
        its corresponding UniSave record as a string. Failed retrievals are
        excluded from the dictionary.
    """
    LOG.info(f'Going to pull {len(entries)} UniSave records')
    
    with logging_redirect_tqdm(loggers = [LOG]):
        unisave_entries = thread_map(pull_from_unisave, entries, 
                                     max_workers = max_workers,
                                     leave = False,
                                     disable = no_progress)
        
    succeeded_unisave_entries = {k:v for k,v in zip(entries, unisave_entries) if v != None}
    
    return succeeded_unisave_entries

