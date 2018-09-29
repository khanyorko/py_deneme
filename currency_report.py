from bs4 import BeautifulSoup as bs
from requests import get

#getting ready to check the list of currencies available on the page
url_base                = "https://kur.doviz.com/serbest-piyasa/"
url_html                = get(url_base)
url_soup                = bs(url_html.content, 'html.parser')

#removing the first element since it has no use for us
base_url                = url_soup.findAll('td', class_ = 'column-row5')[1:]
urls_curr               = list(map(lambda x: x.find('a').attrs['href'].split('/')[-1], url_soup.findAll('td', class_ = 'column-row5')[1:]))

#creating a dict with the same size of currencies we can access
doviz = { x : {
        'url_add' : x,                  # the part we add to the url_base for accessing the currency's page
        'curr_name' : ' ',              # name of the currency
        'source'    : ' ',              # source of the data
        'buy'       : 0,          
        'sell'      : 0,
        'update'    : ' '               # last time data was updated
        } for x in urls_curr}

def getData(url):
    page_url            = url
    page_html           = get(page_url)
    page_soup           = bs(page_html.content, 'html.parser')
    
    curr_name           = page_soup.find('span', class_ = 'left').get_text()
    source              = page_soup.find('span', class_ = 'btn btn-third').get_text()
    soup_data           = list(map(lambda x: x.get_text(), page_soup.find('div', class_ = 'data').findChildren()))
    
    soup_data[2]        = float(soup_data[2].replace(',', '.'))
    soup_data[5]        = float(soup_data[5].replace(',', '.'))
    soup_data[-1]       = soup_data[-1].split(' ')[-1]
    
    return curr_name, source, soup_data[2], soup_data[5], soup_data[-1]
    
def refresh():
    for curr in doviz.keys():
        data = getData(url_base + doviz[curr]['url_add'])
        
        doviz[curr]['curr_name']  = data[0]
        doviz[curr]['source']     = data[1]
        doviz[curr]['buy']        = data[2]
        doviz[curr]['sell']       = data[3]
        doviz[curr]['update']     = data[4]
refresh()
