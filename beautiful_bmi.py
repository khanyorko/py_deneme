from requests import get
from bs4 import BeautifulSoup


def getBMI(age, gender, height, weight):
    url_        = "https://www.calculator.net/bmi-calculator.html?ctype=metric&cage={}&csex={}&cheightmeter={}&ckg={}".format(age, gender, height, weight)
                #  https://www.calculator.net/bmi-calculator.html?ctype=metric&cage=99&csex=f&cheightmeter=66&ckg=99
    page        = get(url_)
    page_soup   = BeautifulSoup(page.content, 'html.parser')
    result      = page_soup.find_all(class_ = 'bigtext')
    return "Your " + result[0].b.get_text().replace('=', 'is'), 'You\'re in ' + result[0].font.get_text() + ' range'
    

# not assuming anything tho
bmi, category = getBMI(99, 'm', 66, 99)
