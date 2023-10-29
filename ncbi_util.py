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

# TaxonInfo object
TaxonInfo = namedtuple('TaxonInfo', ['rank', 'sci_name', 'taxon_id', 'lineage'])

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

def default_preprocessor(raw_name: str):
    """ Input:
            raw_name: str - the name of an organism. i.e. 'Asellus_aquaticus', 
                'Chelifera', 'Ephemerella_aroni_aurivillii', etc.
                The point of this 
        Output:
            str - the name of the organism, processed to obviate the issues I 
                ran into with the dataset I'm working with right now:
                    - organism names that end with "_sp", "_adult", or "_larva"
                    - organism names that contain more than 2 parts (this may be
                        due to ambiguity - I might reexamine this later)
    
    This is an example of a preprocessor for the NCBI class. My particular use
    case was processing directory names. 
    This isn't perfect, nor is it magic.
    It makes assumptions that might be false (that's why I'm calling it an 
    "example preprocessor"). It's moderately likely that you might need to 
    write your own processor that's keyed to the idiosyncracies of the 
    particular dataset you're working with.
    If there are misspelled names you'll have to go in and change them.

    TODO: think about ading some sensible "space between words" handling other
    than underscores ('+','-', ' ', etc.). Maybe camelcase too
    """
    # if there are numbers in the name, we remove them
    raw_name = re.sub('[0-9]', '', raw_name)
    # gets rid of some extra bits on the end
    # NOTE: might want to add a "suffixes" argument or something to make this
    # configurable
    raw_name = re.sub('_sp$|_adult$|_larva$', '', raw_name)
    parts = raw_name.split('_')
    if len(parts) > 1:
        return '+'.join([parts[0], parts[-1]])
    else:
        return parts[0]

class NCBI:
    def __init__(
            self,
            email: str = EMAIL,
            tool: str = TOOL,
            preprocessor: func = None,
            return_ranks: list = None,
            max_attempts: int = 3,
            timeout: int = 10
            ):
        """ Input:
                email: str - the email address that will be used when making the
                    request to the NCBI API. Must be a valid email with no
                    internal spaces.
                tool: str - the name of application making the E-utility call.
                    Value must be a string with no internal spaces.
                preprocessor: func - a function to pre-process the organism
                    names being fed to this. If None, the default preprocessor
                    is used instead
                return_ranks: list - the ranks that will be kept and returned
                max_attempts: int - the number of retries to make  to the API
                    before giving up
                timeout: int - the length of time in seconds to wait before 
                    assuming something went wrong with the request

        Information on parameters, syntax, etc. for the API (including the
        "tool" and "email" parameters for this class) can be found here:
            https://www.ncbi.nlm.nih.gov/books/NBK25499/

        """
        self.email = check_ncbi_param(email, 'email')
        self.tool = check_ncbi_param(tool, 'tool')
        self.return_ranks = return_ranks
        self.disambiguate = []
        self.no_match = []
        self.organisms_known = dict()
        if not preprocessor:
            self.preprocessor = default_preprocessor

    def make_req(
            self,
            url: str,
            payload: dict,
            ):
        """ Input:
                url: str - the url we want to make the request to
                payload: dict - contains the values we want to pass as parameters
                    in the URL's query string
                
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
        return id_list

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
                taxon_info: TaxonInfo namedtuple - the keys are:
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
        return TaxonInfo(**taxon_info)

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
        taxid_list = self.organism_to_id(organism, verbose)

        # If more than 1 taxon id is returned that means we can't match the name
        # given to a single organism. 
        # The user will need to disambiguate the name
        if len(taxid_list) > 1:
            self.disambiguate.append(organism)
            return False

        # If no taxon ids were returned, that means we weren't able to find a
        # match. The user will have to check the spelling of the name.
        if not taxid_list:
            self.no_match.append(organism)
            return False
        
        if not taxid_list[0].isdigit():
            raise ValueError(f'Expected API to return one taxon id consisting \
                            of all decimal characters. Returned {id_list[0]} \
                            instead.')
        taxid = int(taxid_list[0])
        tree = self.etree_from_id(taxid)
        tax_dict = self.etree_to_dict(tree)
        return tax_dict

    def match(
        self,
        org_name: str,
        verbose: bool = False):
        """ Input:
                org_name: str - the name of the organism we want to try to
                    find a match for in the NCBI database.
                verbose: bool - if True, some information will be printed as
                    the data is retrieved.
            Output:
                taxon_info: TaxonInfo - a named tuple containing the organism's:
                        - rank
                        - scientific name
                        - taxonomic ID
                        - lineage 
                    If a match cannot be made, False is returned instead.
        """
        # skips the call to the API if the data for the organism has already
        # been retrieved
        if org_name in self.organisms_known:
            return self.organisms_known[org_name]

        taxon_info = self.organism_to_dict(org_name)
        self.organisms_known[org_name] = taxon_info
        return taxon_info

    def fix(self,
        org_name: str,
        verbose: bool = False):

        pass
