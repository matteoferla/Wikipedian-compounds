# Wikipedian-compounds
Parsing Wikipedia Chembox and Reactionbox data for fun.
Actually, for reactions the view count and the page length might be useful for ranking reactions.

English Wipedia pages with the templates `['drugbox', 'chembox','Drugbox', 'Chembox']`,
and separately with `['Reactionbox', 'reactionbox']` and these are used to tell compounds and reactions from countries.

## Parsing code
Parsing was done as previous projects: see [parsing notes](parsing.md).

SMILES are a bit of a mess...
```python
import re

smiles_pattern = re.compile(r"smiles\s*=\s*(.*?)[\{|\n]", re.IGNORECASE)
def combine_smilestrings(row: pd.Series) -> str:
    for k in ('smiles', 'SMILES', 'SMILES1', 'SMILES2', 'SMILES3', 'SMILES4', 'SMILES5'):
        if isinstance(row[k], str) and row[k].strip():
            return row[k].strip()
    rex = re.search(smiles_pattern, compound_pages[row['title']])
    if rex:
        return rex.group(1)
    return ''
              
compounds['combined_SMILES'] = compounds.apply(combine_smilestrings, axis=1)
```
Elements seem to have a weird infobox for H&S info and these are highly read, so marking elements is required.
```python
elements = ['Actinium', 'Aluminum', 'Americium', 'Antimony', 'Argon', 'Arsenic', 'Astatine', 'Barium', 'Berkelium', 'Beryllium', 'Bismuth', 'Bohrium', 'Boron', 'Bromine', 'Cadmium', 'Calcium', 'Californium', 'Carbon', 'Cerium', 'Cesium', 'Chlorine', 'Chromium', 'Cobalt', 'Copper', 'Curium', 'Darmstadtium', 'Dubnium', 'Dysprosium', 'Einsteinium', 'Erbium', 'Europium', 'Fermium', 'Fluorine', 'Francium', 'Gadolinium', 'Gallium', 'Germanium', 'Gold', 'Hafnium', 'Hassium', 'Helium', 'Holmium', 'Hydrogen', 'Indium', 'Iodine', 'Iridium', 'Iron', 'Krypton', 'Lanthanum', 'Lawrencium', 'Lead', 'Lithium', 'Livermorium', 'Lutetium', 'Magnesium', 'Manganese', 'Meitnerium', 'Mendelevium', 'Mercury', 'Molybdenum', 'Moscovium', 'Neodymium', 'Neon', 'Neptunium', 'Nickel', 'Nihonium', 'Niobium', 'Nitrogen', 'Nobelium', 'Oganesson', 'Osmium', 'Oxygen', 'Palladium', 'Phosphorus', 'Platinum', 'Plutonium', 'Polonium', 'Potassium', 'Praseodymium', 'Promethium', 'Protactinium', 'Radium', 'Radon', 'Rhenium', 'Rhodium', 'Roentgenium', 'Rubidium', 'Ruthenium', 'Rutherfordium', 'Samarium', 'Scandium', 'Seaborgium', 'Selenium', 'Silicon', 'Silver', 'Sodium', 'Strontium', 'Sulfur', 'Tantalum', 'Technetium', 'Tellurium', 'Tennessine', 'Terbium', 'Thallium', 'Thorium', 'Thulium', 'Tin', 'Titanium', 'Tungsten', 'Ununbium', 'Ununhexium', 'Ununoctium', 'Ununpentium', 'Ununquadium', 'Ununseptium', 'Ununtrium', 'Uranium', 'Vanadium', 'Xenon', 'Ytterbium', 'Yttrium', 'Zinc', 'Zirconium']

compounds['is_element'] = compounds['title'].isin(elements)
```
## Top 50 compounds
What are the top 50 most read-about compounds? Monthly views.
```python
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display, HTML

scomps = compounds.loc[compounds.combined_SMILES.astype(bool)]
maximum = 50
svg: str = Draw.MolsToGridImage(scomps.combined_SMILES[:maximum].apply(Chem.MolFromSmiles).to_list(),
                     legends=(scomps.title[:maximum] + ': '+(scomps.counts_2022H2/6)[:maximum].astype(int).astype(str)).to_list(),
                    molsPerRow=5,
                        useSVG=True)
with open('top.svg', 'w') as fh:
    fh.write(svg)
    
display(HTML(svg))
```
![top 50 most read-about compounds](top.svg)

The limonene is a weird 2022 trend: [pageview analytics](https://pageviews.wmcloud.org/?project=en.wikipedia.org&platform=all-access&agent=user&redirects=0&start=2015-07&end=2023-01&pages=Paracetamol%7CLimonene%7CCaffeine%7CEthanol%7CMDMA)

Of the [20,638 compounds](compounds.csv), [930 have had zero readers](unread.txt). Some may be parsing errors, like `Tröger's base` whose umlaut made the matching fail, but they seem mostly to obscure salts.

## Reactions

[345 reactions](reactions.csv). Suzuki is most popular, while Claisen is #44 and Hantzsch is #48? Come on Wikipedia users!