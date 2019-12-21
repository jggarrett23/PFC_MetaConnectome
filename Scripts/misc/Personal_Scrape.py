import bs4 
import requests
from bs4 import BeautifulSoup as soup
import urllib
import pandas as pd
import re
import openpyxl

#Act like actual browser
headers = {'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/53.0.2785.143 Safari/537.36'}

#Input query for search. Later, make search url a customizable input or list of inputs.
query='exercise cognition AND free full text[sb]'
q = urllib.quote_plus(query)
search_url='https://www.ncbi.nlm.nih.gov/pubmed/?term='+ q

#Get search url and parse webpage.
results=requests.get(search_url, headers=headers, timeout=5)
if results.status_code != 200:
    print False
else:
    page = soup(results.content, "lxml")
    results.close()

#Find all 'div' classes that contain the result url and place into a "container" variable.
result_containers=page.findAll("div", {"class": "rprt"})

#Grab url from result url container by searching for 'a','href' tags, and place into list.
results_url_list=[]
for searchWrapper in result_containers:
    article_link=searchWrapper.find('a')["href"].encode('ascii','ignore')
    new_url=search_url.replace(q,'')[:28]+article_link
    results_url_list.append(new_url)

"""Navigate to result url, and grab content from its page (in this case the abstract).
    In upper right hand corner of Pubmed result page should be link to full text. Searches
    for it by looking for 'div' tag of icon, and places link into list."""
    
paper_link=[]
for new_link in results_url_list:
    abstract=requests.get(new_link, headers=headers, timeout=5)
    if abstract.status_code !=200:
        print False
    else:
        abstract_page=soup(abstract.content, "lxml")
        abstract.close()
        full_text_container=abstract_page.findAll("div", {"class": "icons portlet"})
        for wrapper in full_text_container:
            full_text_link=wrapper.find('a')["href"]
            paper_link.append(full_text_link.encode('ascii','ignore'))

#Navigate to full text page.
new_results={}
try:
    for link in paper_link:
        new_results[link]=requests.get(link, headers=headers, timeout=5)
        new_results[link].close()
except requests.exceptions.ConnectionError:
    pass

dic={}
working_links=[]
heading_container={}
table_container={}
for link in new_results.keys():
    if new_results[link].status_code != 200:
        print False, link, new_results[link].status_code
    else:
        #Get full text content.
        full_text_page=soup(new_results[link].content, "lxml")
        #Get the title of the article.
        title=full_text_page.title.text.strip().encode('ascii','ignore')
        #Use regular expression to find all of the headings of the article
        heading_container[link]=full_text_page.findAll(re.compile(r'h\d+'))
        #Use regular expression to find all tables in the article
        table_container[link]=full_text_page.findAll(re.compile(r'td'))
        
        """If word 'Results' in the text of the heading or table tag, place
        into dict with link as key and text of heading as value."""
        
        for x in heading_container[link]:
            if re.search('Results',x.text, re.IGNORECASE):
                dic[title]={link: x.parent.get_text().strip().encode('ascii','ignore')}
                working_links.append(link)
            if 'Results' in x:
                dic[title]={link: x.parent.get_text().strip().encode('ascii','ignore')}
                working_links.append(link)
        for table in table_container[link]:
            if re.search('Results', table.text, re.IGNORECASE):
                dic[title]={link: table.parent.get_text().strip().encode('ascii','ignore')}
                working_links.append(link)
            if 'Results' in table:
                dic[title]={link: table.parent.get_text().strip().encode('ascii','ignore')}
                working_links.append(link)
#Return number of links unable to get a 'Results' paragraph from.
number_of_failures=len(new_results.keys())-len(set(working_links))
print 'Unable to get data from %s links.'%number_of_failures

#Create pandas dataframe to export data to.
dataset=pd.DataFrame(columns=['Title','Link','Text'])
dataset['Title']=dic.keys()
for title in dic.keys():
    for n in range(len(dataset)):
        if title == dataset['Title'][n]:
            dataset.loc[n, 'Link']=dic[title].items()[0][0]
            dataset.loc[n, 'Text']=dic[title].items()[0][1]
            
#Use Openxl to load excel work book with predefined formulas to filter 'Results' paragraph further
book = openpyxl.load_workbook('/Users/owner/Desktop/Article_Database 2.xlsx')
writer = pd.ExcelWriter('/Users/owner/Desktop/Article_Database 2.xlsx', engine='openpyxl') 
writer.book = book
writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
dataset.to_excel(writer, 'Sheet1', index=False)
writer.save()