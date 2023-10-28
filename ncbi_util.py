import os
import time
import re
import json
from pathlib import Path
import requests
# because xml is dangerous
from defusedxml import ElementTree as ET
from xml.etree.ElementTree import Element
import pandas as pd
from collections import namedtuple, defaultdict
from tqdm.notebook import tqdm


BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

EMAIL = os.environ['NCBI_EMAIL_ADDR']
TOOL = os.environ['NCBI_TOOL_NAME']
# The length of time we wait between making requests to the API.
# The docs imply that you can make up to 3 requests per second without an API
# key, but .5 is a bit safer (and more polite)
SLEEP_INTERVAL = .5
RETURN_RANKS = ['order', 'family', 'genus', 'species']
JSON_PATH = 'data/IDA/json'

def check_ncbi_param(ncbi_param: str, param_name: str):
    """ Input:
            ncbi_param: str - the value of a parameter we're going to pass to
                the NCBI API ("email" or "tool").
            param_name: str - the name of the NCBI parameter
        Output:
            Either reurns ncbi_param unaltered, or raises a value error.

    The NCBI API docs specify that the parameters "email" and "tool" should be
    strings with no internal spaces. I think it's probably fine (they get
    encoded by the requests library - so spaces become "+"), maybe I'll get rid
    of this later when I can test more.
    """
    if (not isinstance(ncbi_param, str)) or (' ' in ncbi_param.strip()):
            raise ValueError(f'"{param_name}" parameter must be a string \
                             containing no internal spaces')
    return ncbi_param

class NCBI:
    def __init__(
            self,
            email: str = EMAIL,
            tool: str = TOOL,
            return_ranks: list = None,
            ):
        """ Input:
                email: str - the email address that will be used when making the
                    request to the NCBI API. Must be a valid email with no
                    internal spaces.
                tool: str - the name of application making the E-utility call.
                    Value must be a string with no internal spaces.

        Information on parameters, syntax, etc. for the API (including the
        "tool" and "email" parameters for this class) can be found here:
            https://www.ncbi.nlm.nih.gov/books/NBK25499/

        """
        self.email = check_ncbi_param(email, 'email')
        self.tool = check_ncbi_param(tool, 'tool')
        self.return_ranks = return_ranks

    def make_req(
            self,
            url: str,
            payload: dict,
            max_attempts: int = 3,
            timeout: int = 10):
        """ Input:
                url: str - the url we want to make the request to
                payload: dict - contains the values we want to pass as parameters
                    in the URL's query string
                max_attempts: int - the number of retries to make before giving up
                timeout: int - the length of time in seconds to wait before assuming
                    something went wrong with the request
            Output:
                req: requests.Response - the response returned by the NCBI API

        This attempts a GET request, and then waits for a bit and tries again if
        it receives a RequestException. It helps if you've got some transient
        network issues between you and the other fella.
        Also: having all the requests go through this method makes it easier
        to keep track of stuff like the slight delay between requests (to avoid
        being rate limited)
        """
        for i in range(max_attempts + 1):
            try:
                req = requests.get(
                    url,
                    params=payload,
                    timeout=timeout)
                req.raise_for_status()
                time.sleep(SLEEP_INTERVAL)
                return req
            except requests.RequestException as e:
                if i == max_attempts:
                    raise e
                print('connection issue. Waiting 1 second and trying again...')
                time.sleep(1)
        raise requests.HTTPError

    def esearch_req(
            self,
            organism: str):
        """ Input:
                organism: str - the name of the organism whose ID we want.
            Output:
                req: requests.Response - the response returned by the esearch
                    endpoint of the NCBI API for the organism in question.

        The NCBI eutils esearch endpoint can return a list of UIDs that match a
        query.
        Docs are here:
        https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
        """
        url = BASE_URL + 'esearch.fcgi'
        payload = {
            'mail': self.email,
            'tool': self.tool,
            'db':'taxonomy',
            'term':organism,
            'rettype':'uilist',
            'retmode':'json'
            }
        req = self.make_req(url, payload)
        return req

    def efetch_req(
            self,
            taxid: int):
        """ Input:
                taxid: int - the taxon ID that we want more information on
            Output:
                req: requests.Response - the response returned by the efetch
                    endpoint of the NCBI API for the taxon in question.

        The NCBI eutils efetch endpoint returns data records for a UID or list of
        UIDs.
        Docs are here:
        https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
        """
        url = BASE_URL + 'efetch.fcgi'
        payload = {
            'mail': self.email,
            'tool': self.tool,
            'db':'taxonomy',
            'id':taxid}
        req = self.make_req(url, payload)
        return req

    def organism_to_id(
            self,
            organism: str,
            return_raw: bool = False,
            verbose: bool = False):
        """ Input:
                organism: str - the name of the organism whose ID we want.
                return_raw: bool - if true, returns the raw idlist returned
                    by the API.
                verbose: bool - if true, prints each organism name before
                    starting the process
            Output:
                int - the taxon id for the given organism
        """
        if verbose:
            print('Starting organism:', organism)
        req = self.esearch_req(organism)
        id_list = req.json()['esearchresult']['idlist']
        if return_raw:
            return id_list
        if len(id_list) > 1:
            raise ValueError(f'Expected API to return one id, but it returned \
                             {len(id_list)} ids instead')
        if not id_list:
            raise ValueError(f'Request returned empty id_list: {req.text}')
        if not id_list[0].isdigit():
            raise ValueError(f'Expected API to return one taxon id consisting \
                             of all decimal characters. Returned {id_list[0]} \
                                instead.')
        return int(id_list[0])

    def etree_from_id(
            self,
            taxid: int):
        """ Input:
                taxid: int - an NCBI taxonomic id.
            Output:
                tree: defusedxml.ElementTree - an element tree containing the
                    lineage information returned by the NCBI efetch API for the
                    organism with the given taxid.
        """
        req = self.efetch_req(taxid)
        tree = ET.fromstring(req.content)
        return tree

    def etree_to_dict(self,
                      root: ET):
        """ Input:
                root: defusedxml.ElementTree - an xml document as returned by 
                    the NCBI API efetch endpoint.
            Output:
                taxon_info: dict - the keys are:
                    - rank: the rank of the organism, i.e. order, phylum, etc.
                    - sci_name: the scientific name of the organism
                    - taxon_id: the taxonomic id (an int)
                    - lineage: a list of dicts. The keys are:
                        - rank:
                        - sci_name: 
                        - taxon_id: 
        Note: Within taxon_info we have the lineage dictionary, which is a list
        of dicts. We're using a list of dicts because some ranks (specifically
        "clade") can appear multiple times.
        I decided to preserve all of the lineage data at this stage and filter
        out unused entries further down the line, so it will be easy to rework -
        just in case I find a use for the extra data somewhere.
        """
        # This gets us the rank, scientific name, and taxon id for the organism
        root_taxon = root.find('Taxon')
        rank, taxon_info = self.parse_taxon_element(root_taxon)
        taxon_info['rank'] = rank
        
        # This gets a list of the taxa in the organism's lineage and then 
        # creates a dictionary of lists, in which each entry is a dict 
        # containing scientific name and taxonomic id
        lineage = defaultdict(list)
        taxa = root_taxon.find('LineageEx').findall('Taxon')
        for taxon in taxa:
            rank, info = self.parse_taxon_element(taxon)
            lineage[rank].append(info)
            
        taxon_info['lineage'] = lineage
        return taxon_info

    def parse_taxon_element(
            self,
            taxon: Element):
        ''' Input:
                taxon: Element - the 'Taxon' element from the element tree which
                    we retrieved from the NCBI API
            Output:
                (rank, info): a tuple. "rank" is the rank of the element (i.e.
                    "order", "phylum", etc.), and "info" is a dictionary
                    containing the scientific name and taxon id for that rank.
        '''
        rank = taxon.find('Rank').text
        info = {
            'sci_name': taxon.find('ScientificName').text,
            'taxon_id': int(taxon.find('TaxId').text)
        }
        return rank, info

    def organism_to_dict(
            self,
            organism: str,
            verbose: bool = False):
        ''' Input:
                organism: str - the organism we're interested in
                verbose: bool - if true, some progress info will be printed as
                    the data is retrieved.
            Output:
                tax_dict: dict - a dictionary containing the taxonomy data about
                    "organism" which was returned by the NCBI API
        '''
        taxid = self.organism_to_id(organism, verbose)
        tree = self.etree_from_id(taxid)
        tax_dict = self.etree_to_dict(tree)
        return tax_dict

    ###
    ### TODO: REMOVE - this should be in a subclass
    ###
    # def preprocess_name(
    #         self,
    #         dir_name: str):
    #     """ Input:
    #             dir_name: str - the name of a directory. This name should correspond
    #                 to the name of an organism. i.e. 'Asellus_aquaticus', 'Chelifera',
    #                 'Ephemerella_aroni_aurivillii', etc.
    #         Output:
    #             str - the name of the organism, processed to obviate the issues I ran
    #                 into with the single dataset I'm working with right now:
    #                     - organism names that end with "_sp", "_adult", or "_larva"
    #                     - organism names that contain more than 2 parts (this may be due
    #                     to ambiguity - I might reexamine this later)
        
    #     This isn't perfect, nor is it magic.
    #     If there are misspelled names you'll have to go in and change them.
    #     """
    #     parts = re.sub('_sp|_adult|_larva', '', dir_name).split('_')
    #     if len(parts) > 1:
    #         return '+'.join([parts[0], parts[-1]])
    #     else:
    #         return parts[0]

    # def get_names_from_dataset(path: str):
    #     """ Input:
    #             path: str - the path to the directory containing the dataset.
    #                 This assumes that the images are stored in directories whose
    #                 names are formatted as "Genus_species".
    #         Output:
    #             organism_list: list of strings - the names of the directories.
    #     """
    #     p = Path(path)
    #     organism_list = [x.name for x in p.iterdir() if x.is_dir()]
    #     return organism_list

    # def dict_from_path(file_path: str):
    #     """ Input:
    #             file_path: str - the path to a file containing json data. 
    #         Output:
    #             a dictionary containing the data in the file specified by 
    #             "file_path". If the file does not exist, an empty dictionary is
    #             returned instead.
    #     """
    #     path = Path(file_path)
    #     if path.exists():
    #         with path.open() as f:
    #             return json.load(f)
    #     else:
    #         return dict()
        
    # def get_taxon_data(taxon_path:str, 
    #                 dir_to_taxon:str, 
    #                 name_to_taxon:str, 
    #                 failed_path:str, 
    #                 preprocessor):
    #     """ Input:
    #             taxon_path: str - the path of the json file in which we want to store
    #                 our taxon data. If the file doesn't exist it will be created.
    #             preprocessor: function - a function that will take a raw directory name
    #                 and handle any weird stuff in it (removing numbers, for exampe)
    #             typo_path: str - the path to a json file containing a dictionary
    #                 that maps the names of directories that have been misspelled to
    #                 corrected versions.
    #                 I'm doing this as a dictionary in a separate json file instead
    #                 of the simpler option of just changing the name of the directory
    #                 manually because I might need to wipe the data - this happened 
    #                 once when an image got corrupted by Tensorflow and it was easier
    #                 to just wipe the directory and unzip it again.
    #                 It's not just a dict at the top of the other json doc for
    #                 similar reasons.
    #         Output:
    #             taxon_data: dict - a dictionary in which the keys are the names of
    #                 directories, and the values are dictionaries containing the
    #                 taxonomic lineage information returned by the API.
    #     """
    #     taxon_data = dict_from_path(taxon_path)
    #     dir2taxon = dict_from_path(dir_to_taxon)
    #     name2taxon = dict_from_path(name_to_taxon)
    #     failed_to_retrieve = list(dict_from_path(failed_path))
        
    #     directory_names = set(get_names_from_dataset(DATA_PATH))
    #     new_directories = directory_names.difference(dir2taxon.keys())
    #     directory_names = list(directory_names)
    #     failed_to_retrieve = list()
    #     print(f'Found {len(directory_names)} directories:')
    #     print(f'\t- {len(directory_names)-len(new_directories)} directories already\
    #     saved to json')
    #     print(f'\t- {len(new_directories)} new directories to retrieve data for\n')
        
    #     progress_bar = tqdm(total=len(new_directories))
    #     for v, dir_name in enumerate(new_directories):
            
    #         organism_name = preprocessor(dir_name)
    #         organism_name = preprocess_name(organism_name)
    # #         print(f'organism_name: {organism_name}')
    #         if organism_name in name2taxon:
    #             dir2taxon[dir_name] = name2taxon[organism_name]
    #             continue
    #         try:
    #             # We want a dictionary where the keys are the names of directories,
    #             # so we only clean up the directory name with preprocess_name() and
    #             # get "organism_name" to make the API call
                
    #             taxon_dict = organism_to_dict(organism_name)
    #             taxon_id = taxon_dict['taxon_id']
    #             name2taxon[organism_name] = taxon_id
    #             dir2taxon[dir_name] = taxon_id
    #             taxon_data[taxon_id] = taxon_dict
    #             progress_bar.reset()
    #             progress_bar.update(v)
    #             progress_bar.set_description(desc=dir_name.rjust(20, '-')[:20])
    #             # Speed isn't important (we're putting in a pause in the function that makes
    #             # the requests anyway, to be considerate to the API we're hitting), so there's
    #             # no reason not to just write to disk on every loop 
    #             with open(taxon_path, 'w+') as f:
    #                 json.dump(taxon_data, f)
    #             with open(dir_to_taxon, 'w+') as f:
    #                 json.dump(dir2taxon, f)
    #             with open(name_to_taxon, 'w+') as f:
    #                 json.dump(name2taxon, f)
    #         except ValueError:
                
    #             progress_bar.reset()
    #             progress_bar.update(v)
    #             progress_bar.set_description(desc=dir_name.rjust(20, '-')[:20])
    #             if organism_name in failed_to_retrieve:
    #                 continue
    #             print(f'Failed to retrieve data for organism "{organism_name}"')
    #             failed_to_retrieve.append(organism_name)
    #             with open(failed_path,'w+') as f:
    #                 json.dump(failed_to_retrieve, f)
    # #         if v > 100:
    # #             break
    #     progress_bar.close()
    #     print(f'\nData for {len(new_directories) - len(failed_to_retrieve)} \
    #         directories successfully retrieved')
    #     print(f'Failed to retrieve {len(failed_to_retrieve)} directories:')
    #     print(failed_to_retrieve)
    #     return taxon_data