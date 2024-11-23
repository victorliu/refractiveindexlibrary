import os
import yaml
import numpy as np
import scipy.interpolate
from io import open

try:
	from yaml import CBaseLoader as BaseLoader
except ImportError:
	from yaml import BaseLoader

class Library:
	def __init__(self, dbpath=None):
		self._index = {}
		if dbpath is not None:
			self.add_database(dbpath)

	def add_database(self, dbpath):
		dbpath = os.path.normpath(dbpath)
		filename = os.path.join(dbpath, 'catalog-nk.yml')
		with open(filename, 'rt', encoding='utf-8') as f:
			catalog = yaml.load(f, Loader=BaseLoader)

		for sh in catalog:
			curshelf = {}
			self._index[sh['SHELF']] = curshelf
			for b in sh['content']:
				if 'DIVIDER' in b: continue
				curbook = {}
				curshelf[b['BOOK']] = curbook
				for p in b['content']:
					if 'DIVIDER' in p: continue
					curbook[p['PAGE']] = {
						'name': p['name'],
						'filename': os.path.join(dbpath, 'data-nk', os.path.normpath(p['data']))
					}

	def get_shelves(self) -> list[str]:
		return list(self._index.keys())

	def get_books(self, shelf: str) -> list[str]:
		if shelf in self._index:
			return list(self._index[shelf].keys())
		return list()

	def get_pages(self, shelf: str, book: str) -> list[str]:
		if shelf in self._index:
			if book in self._index[shelf]:
				return list(self._index[shelf][book].keys())
		return list()

	def get_entry(self, shelf: str, book: str, page: str) -> dict:
		if shelf in self._index:
			if book in self._index[shelf]:
				if page in self._index[shelf][book]:
					return self._index[shelf][book][page]
		return None

	def get_material(self, shelf: str, book: str, page: str):
		entry = self.get_entry(shelf, book, page)
		if entry is not None:
			return Material(entry['filename'])
		else:
			return None



class Material:
	""" Material class"""

	def __init__(self, filename):
		"""

		:param filename:
		"""
		self.refractiveIndex = None
		self.extinctionCoefficient = None
		self.originalData = None

		with open(filename, "rt", encoding="utf-8") as f:
			material = yaml.load(f, Loader=BaseLoader)

		for data in material['DATA']:
			if (data['type'].split())[0] == 'tabulated':
				rows = data['data'].split('\n')
				splitrows = [c.split() for c in rows]
				wavelengths = []
				n = []
				k = []
				for s in splitrows:
					if len(s) > 0:
						wavelengths.append(float(s[0]))
						n.append(float(s[1]))
						if len(s) > 2:
							k.append(float(s[2]))

				if (data['type'].split())[1] == 'n':
					if self.refractiveIndex is not None:
						Exception('Bad Material YAML File')

					self.refractiveIndex = RefractiveIndexData.setupRefractiveIndex(
						formula=-1,
						wavelengths=wavelengths,
						values=n
					)
					self.originalData = {
						'wavelength (um)': np.array(wavelengths), 
						'n' : np.array(n)
					}
				elif (data['type'].split())[1] == 'k':
					self.extinctionCoefficient = ExtinctionCoefficientData.setupExtinctionCoefficient(wavelengths, n)
					self.originalData = {
						'wavelength (um)': np.array(wavelengths), 
						'n' : 1j*np.array(n)
					}

				elif (data['type'].split())[1] == 'nk':
					if self.refractiveIndex is not None:
						Exception('Bad Material YAML File')

					self.refractiveIndex = RefractiveIndexData.setupRefractiveIndex(
						formula=-1,
						wavelengths=wavelengths,
						values=n
					)
					self.extinctionCoefficient = ExtinctionCoefficientData.setupExtinctionCoefficient(wavelengths, k)
					self.originalData = {
						'wavelength (um)': np.array(wavelengths), 
						'n' : np.array(n) + 1j*np.array(k)
					}

			elif (data['type'].split())[0] == 'formula':

				if self.refractiveIndex is not None:
					Exception('Bad Material YAML File')

				formula = int((data['type'].split())[1])
				coefficents = [float(s) for s in data['coefficients'].split()]
				for k in ['range','wavelength_range']:
					if k in data:
						break
				rangeMin = float(data[k].split()[0])
				rangeMax = float(data[k].split()[1])

				self.refractiveIndex = RefractiveIndexData.setupRefractiveIndex(
					formula=formula,
					rangeMin=rangeMin,
					rangeMax=rangeMax,
					coefficients=coefficents
				)
				wavelengths = np.linspace(rangeMin, rangeMax, 1000)
				self.originalData = {
					'wavelength (um)': wavelengths,
					'n' : self.refractiveIndex.getRefractiveIndex(wavelength_um=wavelengths)
				}
				

	def get_n(self, wavelength_um, bounds_error=True):
		"""

		:param wavelength:
		:return: :raise Exception:
		"""
		if self.refractiveIndex is None:
			raise Exception('No refractive index specified for this material')
		else:
			return self.refractiveIndex.getRefractiveIndex(wavelength_um, bounds_error=bounds_error)

	def get_k(self, wavelength_um, bounds_error=True):
		"""

		:param wavelength:
		:return: :raise NoExtinctionCoefficient:
		"""
		if self.extinctionCoefficient is None:
			raise NoExtinctionCoefficient('No extinction coefficient specified for this material')
		else:
			return self.extinctionCoefficient.getExtinctionCoefficient(wavelength_um, bounds_error=bounds_error)

	def get_nk(self, wavelength_um, bounds_error=True):
		n = 0
		k = 0
		if self.refractiveIndex is not None:
			n = self.refractiveIndex.getRefractiveIndex(wavelength_um, bounds_error=bounds_error)
		if self.extinctionCoefficient is not None:
			k = self.extinctionCoefficient.getExtinctionCoefficient(wavelength_um, bounds_error=bounds_error)
		return n + 1j * k

	def get_wavelength_range(self) -> tuple[float, float]:
		return (self.refractiveIndex.rangeMin, self.refractiveIndex.rangeMax)

class RefractiveIndexData:
	"""Abstract RefractiveIndex class"""

	@staticmethod
	def setupRefractiveIndex(formula, **kwargs):
		"""

		:param formula:
		:param kwargs:
		:return: :raise Exception:
		"""
		if formula >= 0:
			return FormulaRefractiveIndexData(formula, **kwargs)
		elif formula == -1:
			return TabulatedRefractiveIndexData(**kwargs)
		else:
			raise Exception('Bad RefractiveIndex data type')

	def getRefractiveIndex(self, wavelength_um):
		"""

		:param wavelength:
		:raise NotImplementedError:
		"""
		raise NotImplementedError('Different for functionally and experimentally defined materials')


class FormulaRefractiveIndexData:
	"""Formula RefractiveIndex class"""

	def __init__(self, formula, rangeMin, rangeMax, coefficients):
		"""

		:param formula:
		:param rangeMin:
		:param rangeMax:
		:param coefficients:
		"""
		self.formula = formula
		self.rangeMin = rangeMin
		self.rangeMax = rangeMax
		self.coefficients = coefficients

	def getRefractiveIndex(self, wavelength_um, bounds_error=True):
		"""

		:param wavelength:
		:return: :raise Exception:
		"""
		wavelength = np.copy(wavelength_um)
		if self.rangeMin <= np.min(wavelength) <= self.rangeMax and self.rangeMin <= np.max(wavelength) <= self.rangeMax or not bounds_error:
			formula_type = self.formula
			coefficients = self.coefficients
			n = 0
			if formula_type == 1:  # Sellmeier
				nsq = 1 + coefficients[0]
				g = lambda c1, c2, w: c1 * (w ** 2) / (w ** 2 - c2 ** 2)
				for i in range(1, len(coefficients), 2):
					nsq += g(coefficients[i], coefficients[i + 1], wavelength)
				n = np.sqrt(nsq)
			elif formula_type == 2:  # Sellmeier-2
				nsq = 1 + coefficients[0]
				g = lambda c1, c2, w: c1 * (w ** 2) / (w ** 2 - c2)
				for i in range(1, len(coefficients), 2):
					nsq += g(coefficients[i], coefficients[i + 1], wavelength)
				n = np.sqrt(nsq)
			elif formula_type == 3:  # Polynomal
				g = lambda c1, c2, w: c1 * w ** c2
				nsq = coefficients[0]
				for i in range(1, len(coefficients), 2):
					nsq += g(coefficients[i], coefficients[i + 1], wavelength)
				n = np.sqrt(nsq)
			elif formula_type == 4:  # RefractiveIndex.INFO
				g1 = lambda c1, c2, c3, c4, w: c1 * w**c2 / (w**2 - c3**c4)
				g2 = lambda c1, c2, w: c1 * w**c2
				nsq = coefficients[0]
				for i in range(1, min(8, len(coefficients)), 4):
					nsq += g1(coefficients[i], coefficients[i+1], coefficients[i+2], coefficients[i+3], wavelength)
				if len(coefficients) > 9:
					for i in range(9, len(coefficients), 2):
						nsq += g2(coefficients[i], coefficients[i+1], wavelength)
				n = np.sqrt(nsq)
			elif formula_type == 5:  # Cauchy
				g = lambda c1, c2, w: c1 * w ** c2
				n = coefficients[0]
				for i in range(1, len(coefficients), 2):
					n += g(coefficients[i], coefficients[i + 1], wavelength)
			elif formula_type == 6:  # Gasses
				n = 1 + coefficients[0]
				g = lambda c1, c2, w: c1 / (c2 - w ** (-2))
				for i in range(1, len(coefficients), 2):
					n += g(coefficients[i], coefficients[i + 1], wavelength)
			elif formula_type == 7:  # Herzberger
				g1 = lambda c1, w, p: c1 / (w**2 - 0.028)**p
				g2 = lambda c1, w, p: c1 * w**p
				n = coefficients[0]
				n += g1(coefficients[1], wavelength, 1)
				n += g1(coefficients[2], wavelength, 2)
				for i in range(3, len(coefficients)):
					n += g2(coefficients[i], wavelength, 2*(i-2))
			elif formula_type == 8:  # Retro
				raise FormulaNotImplemented('Retro formula not yet implemented')
			elif formula_type == 9:  # Exotic
				raise FormulaNotImplemented('Exotic formula not yet implemented')
			else:
				raise Exception('Bad formula type')

			n = np.where((self.rangeMin<=wavelength) & (wavelength<=self.rangeMax), n, np.nan)
			return n
		else:
			raise Exception('Wavelength {} um is out of bounds. Correct range: ({} um, {} um)'.format(wavelength, self.rangeMin, self.rangeMax))


class TabulatedRefractiveIndexData:
	"""Tabulated RefractiveIndex class"""

	def __init__(self, wavelengths, values):
		"""

		:param wavelengths:
		:param values:
		"""
		self.rangeMin = np.min(wavelengths)
		self.rangeMax = np.max(wavelengths)

		if self.rangeMin == self.rangeMax:
			self.refractiveFunction = values[0]
		else:
			self.refractiveFunction = scipy.interpolate.interp1d(wavelengths, values, bounds_error=False)

	def getRefractiveIndex(self, wavelength_um, bounds_error=True):
		"""

		:param wavelength:
		:return: :raise Exception:
		"""
		wavelength = np.copy(wavelength_um)
		if self.rangeMin == self.rangeMax and self.rangeMin == wavelength:
			return self.refractiveFunction
		elif self.rangeMin <= np.min(wavelength) <= self.rangeMax and self.rangeMin <= np.max(wavelength) <= self.rangeMax and self.rangeMin != self.rangeMax or not bounds_error:
			return self.refractiveFunction(wavelength)
		else:
			raise Exception('Wavelength {} um is out of bounds. Correct range: ({} um, {} um)'.format(wavelength, self.rangeMin, self.rangeMax))


class ExtinctionCoefficientData:
	"""ExtinctionCofficient class"""

	@staticmethod
	def setupExtinctionCoefficient(wavelengths, values):
		"""

		:param wavelengths:
		:param values:
		:return:
		"""
		return ExtinctionCoefficientData(wavelengths, values)

	def __init__(self, wavelengths, coefficients):
		"""

		:param wavelengths:
		:param coefficients:
		"""
		self.extCoeffFunction = scipy.interpolate.interp1d(wavelengths, coefficients, bounds_error=False)
		self.rangeMin = np.min(wavelengths)
		self.rangeMax = np.max(wavelengths)

	def getExtinctionCoefficient(self, wavelength_um, bounds_error=True):
		"""

		:param wavelength:
		:return: :raise Exception:
		"""
		wavelength = np.copy(wavelength_um)
		if self.rangeMin <= np.min(wavelength) <= self.rangeMax and self.rangeMin <= np.max(wavelength) <= self.rangeMax or not bounds_error:
			return self.extCoeffFunction(wavelength)
		else:
			raise Exception('Wavelength {} um is out of bounds. Correct range: ({} um, {} um)'.format(wavelength, self.rangeMin, self.rangeMax))


class FormulaNotImplemented(Exception):
	def __init__(self, value):
		self.value = value

	def __str__(self):
		return repr(self.value)


class NoExtinctionCoefficient(Exception):
	def __init__(self, value):
		self.value = value

	def __str__(self):
		return repr(self.value)
	
