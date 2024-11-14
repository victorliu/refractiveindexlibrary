# Python interface to RefractiveIndex database

This package provides a programmatic interface to the [refractiveindex.info](https://refractiveindex.info) database, as well as allowing users to import their own libraries in the same database format.

This package is built off of the [refractiveindex](https://github.com/toftul/refractiveindex) package.

## Installation

```
pip install refractiveindexlibrary
```

In order to use the refractiveindex.info data, you must also clone or download the [database repository](https://github.com/polyanskiy/refractiveindex.info-database).

## Usage


```python
from refractiveindexlibrary import Library

lib = Library()

# Import the refractiveindex.info database
lib.add_database(dbpath = os.path.join(os.path.expanduser('~'), 'refractiveindex.info-database', 'database'))

shelves = lib.get_shelves()
print(shelves)

shelf = shelves[0]
books = lib.get_books(shelf)
print(books)

book = books[0]
pages = lib.get_pages(shelf, book)
print(pages)


SiO = lib.get_material(shelf='main', book='SiO', page='Hass')

wavelength_um = 0.600

SiO.get_n(wavelength_um)
# (1.96553846)

SiO.get_k(wavelength_um)
# (0.001)
```
