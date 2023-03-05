```python
import mwxml, bz2
from wiki_category_analyser import WikicatParser

dump = mwxml.Dump.from_file(bz2.open('dumps/enwiki-latest-pages-articles.xml.bz2', mode='rt', encoding='utf-8'))

print(dump.site_info.name, dump.site_info.dbname)

# ----------------

import re, functools

data = []


compound_pages = {}
for page in dump:
    revision = next(page)
    if revision.text is None:
        continue
    for k in ('drugbox', 'chembox','Drugbox', 'Chembox'):
        if k in revision.text:
            break
    else:
        continue
    compound_pages[page.title] = revision.text
    
# ----------------
    
import gzip, pickle

with gzip.open('compound_pages.pkl.gz', 'wb') as fh:
    pickle.dump(compound_pages, fh)
    
# ----------------

import gzip
import pandas as pd


wp = WikicatParser('', wanted_templates=Engli)

data = []
for title, page in compound_pages.items():
    data.append({'title': title, **wp.parse_templates(page)})
    
compounds = pd.DataFrame(data)

with gzip.open('compounds.pkl.gz', 'wb') as fh:
    compounds.to_pickle(fh)
  
# ----------------

# this works only on my setup sorry...
import json, gzip
for date_tag in ('202208', '202209',  '202210', '202211', '202212', '202301'):
    with gzip.open(f'extracted_dumps/pageviews-{date_tag}.json.gz', 'rb') as fh:
        counts = json.load(fh)
        compounds[f'counts_{date_tag}'] = compounds.title.str.replace(' ', '_').apply(counts.get).fillna(0).astype(int)
        
# second half of 2022...
compounds[f'counts_2022H2'] = 0

for date_tag in ('202208', '202209',  '202210', '202211', '202212', '202301'):
    compounds[f'counts_2022H2'] += compounds[f'counts_{date_tag}']
    
compounds = compounds.sort_values('counts_2022H2', ascending=False)
```