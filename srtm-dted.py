# This software is public domain.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# Written by Philipp Lay <philipp.lay@illunis.net>.
#

from math import trunc
from struct import unpack


class ParseError(Exception): pass



# decode 16bit signed big-endian non-complemented as used
# as elevation data in DTED files.
# if x == 0xffff it is, acording to our spec, interpreted as invalid/unknown,
# and 'None' is returned instead of -32767. This should prevent using this
# number by accident in a calculation.
#
def decode_be16nc(x):
	if x == 0xffff:
		return None
	elif x >= 0x8000:
		return 0x8000-x

	return x


class GeoCoord:
	type = None

	# hemi is either 'N' or 'S' for type=='latitude'
	# or 'E' / 'W' for type == 'longitude'
	hemi = None

	# numerical coordinate in hms
	hms = None


	def __init__(self, str):
		self.from_str(str)

	# format: D.MMSSH
	# where H is the hemisphere
	# for example: 47.2516N
	# for 47 deg 25 min and 16 sec northern hemisphere
	def from_str(self, str):
		H = str[-1].upper()

		if H == 'N' or H == 'S':
			self.hemi = H
			self.type = 'latitude'
		elif H == 'W' or H == 'E':
			self.hemi = H
			self.type = 'longitude'


		self.hms = float(str[0:-1])


	def __str__(self):
		return str(self.hms) + self.hemi

	def deci(self):
		x = self.hms
		erg = trunc(x)

		for i in range(2):
			erg *= 60
			x *= 100
			erg += trunc(x) % 100

		return erg / 3600


	def sub_coord(self):
		x0 = (trunc(100 * self.hms) % 100) % 15
		x1 = trunc(10000 * self.hms) % 100

		if self.hemi == 'N' or self.hemi == 'E':
			return 60*x0 + x1
		else:
			return 900 - 60*x0 - x1



	# give the name (-component) of the map this coordinate belongs to
	def mapname(self):
		x = trunc(self.hms * 100)

		if self.hemi == 'N' or self.hemi == 'E':
			x = ( x - ((x % 100) % 15) ) * 100
		else:
			x = ( x - ((x % 100) % 15) + 15 ) * 100

		if self.type == 'latitude':
			return self.hemi + '%06d' % x
		elif self.type == 'longitude':
			return self.hemi + '%07d' % x




class DTED:
	data = None
	mapname = None

	def __init__(self, file):
		self.data = []
		self.fromfile(file)


	def fromfile(self, file):
		self.mapname = file
		with open(file, 'rb') as f:
			#
			# parse UHL header
			#

			# check magic and version
			if f.read(3) != b'UHL':
				raise ParseError('wrong magic')
			if f.read(1) != b'1':
				raise ParseError('wrong version')

			self.longitude = f.read(8).decode('utf-8')
			self.latitude = f.read(8).decode('utf-8')

			self.longival = int(f.read(4))/10
			self.latival = int(f.read(4))/10

			temp = f.read(4)
			if temp == b'NA  ':
				self.vacc = None
			else:
				self.vacc = int(temp)

			self.usc = f.read(3)
			self.uref = f.read(12)

			self.num_long = int(f.read(4))
			self.num_lat = int(f.read(4))

			self.mult_acc = int(f.read(1))

			self.reserved = f.read(24)


			#
			# skip other headers
			#
			f.seek(648, 1)		# DSI
			f.seek(2700, 1)		# ACC


			#
			# now go for the data blocks
			#

			for i in range(self.num_long):
				# check magic number of record
				if int(unpack('>B', f.read(1))[0]) != 0xaa:
					raise ParseError('wrong magic in data block')

				seq = unpack('>I', b'\x00' + f.read(3))[0]

				long_cnt = unpack('>H', f.read(2))[0]
				if long_cnt != i:
					raise ParseError('unexpected longitude number: ' + str(long_cnt))

				lat_cnt = unpack('>H', f.read(2))[0]
				if lat_cnt != 0:
					raise ParseError('latitude count not zero')

				# read elevations
				rowdata = f.read(2 * self.num_lat)

				# check values with checksum
				# (checksum is calculated in regular 2-complement)
				checksum = unpack('>i', f.read(4))[0]
				rowsum = sum(unpack('>' + 901*'h', rowdata))
				if rowsum != checksum:
					raise ParseError('checksum failed on longitude ' + str(i) + '., should be: ' + str(checksum) + ' but is ' + str(rowsum))

				# add row to matrix, this time values are
				# interpreted as signed big endian 16it non-complement
				self.data.append([ decode_be16nc(x) for x in unpack('>' + self.num_lat * 'H', rowdata) ])



def height(coordstr):
	c0, c1 = ( GeoCoord(s) for s in coordstr.split(' ') )

	if c0.type == 'latitude' and c1.type == 'longitude':
		lat, long = c0, c1
	elif c0.type == 'longitude' and c1.type == 'latitude':
		lat, long = c1, c0
	else:
		raise Exception('you need longitude and latitude')

	dted = DTED(long.mapname() + lat.mapname() + '_SRTM_1_DEM.dt2')

	return dted.data[long.sub_coord()][lat.sub_coord()]



# Example:
#
# Print height of the german Zugspitze, which is
# at '47.2516N 10.5911E'
#
# you will need the corresponding DTED map for this:
#   E0104500N471500_SRTM_1_DEM.dt2
#

#coord = '47.2516N 10.5911E';
#print(coord, " => ", height(coord));
