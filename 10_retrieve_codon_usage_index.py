from bs4 import BeautifulSoup
import requests
import utils

def main():
    for host, url in utils.codon_usage_index_urls.items(): # Loops through dictionary containing host name and URL
        write_path = "data/codon_usage_database/{0}.txt".format(host) # Destination of output fil
        condon_usage_index = get_codon_usage_index(url) # Codon usage index scraped from web data using URL
        utils.write_txt(write_path, condon_usage_index)
        print("Retrieved {0} codon usage index.".format(host))

# Scrapes web data from codon usage database to return codon usage index
def get_codon_usage_index(url):
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'html.parser')
    return soup.text
    
if __name__ =="__main__":
    main()
