from bs4 import BeautifulSoup as bs
from requests import get

def getSoup(url):
    page_url    = url
    page_cont   = get(page_url)
    return bs(page_cont.content, 'html.parser')

def getData(soup):
    raw_data    = list(map(lambda x: x.find_all('td'), soup))
    dict_data   =   {
                        day[0].get_text().split('2018')[0] + '2018' :
                        {
                          'weather'     : day[2].get_text(),
                          'day'         : day[3].get_text(),
                          'night'       : day[4].get_text()
                        }
                        for day in raw_data
                    }
                    
    return dict_data
    

base_url        = "https://havadurumu.com.tr/turkiye-sehir-listesi/"
base_soup       = getSoup(base_url)

city_list       = base_soup.find('ul', class_ = 'list-unstyled city_list').find_all('a')
city_links      =   [ 
                        [x.get_text(), "https://havadurumu.com.tr/havadurumu/" + x['href'].split('/')[-1]]
                        for x in city_list
                    ]
#serving 81 cities of Turkey in as a dict.                    
city_dict       =   {
                        city[0] : getData(getSoup(city[1]).find('table', class_ = 'days15').find_all('tr')[1:])
                        for city in city_links
                    }
