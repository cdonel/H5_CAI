from bs4 import BeautifulSoup
import requests
import utils

def run(_host=None):
    print('Running: Retrieving codon usage index.')
    if _host == None:
        for host, url in utils.codon_usage_index_urls.items(): # Loops through dictionary containing host name and URL
            write_path = "data/codon_usage_database/{0}.txt".format(host) # Destination of output fil
            condon_usage_index = request_codon_usage_index(url) # Codon usage index scraped from web data using URL
            utils.write_txt(write_path, condon_usage_index)
            print("Retrieved {0} codon usage index.".format(host))
    else:
        for host, url in utils.codon_usage_index_urls.items(): # Loops through dictionary containing host name and URL
            if host == _host:
                write_path = "data/codon_usage_database/{0}.txt".format(host) # Destination of output fil
                condon_usage_index = request_codon_usage_index(url) # Codon usage index scraped from web data using URL
                utils.write_txt(write_path, condon_usage_index)
                print("Retrieved {0} codon usage index.".format(host))

# Scrapes web data from codon usage database to return codon usage index
def request_codon_usage_index(url):
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'html.parser')
    return soup.text
    
