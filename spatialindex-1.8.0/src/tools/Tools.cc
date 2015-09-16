/******************************************************************************
 * Project:  libspatialindex - A C++ library for spatial indexing
 * Author:   Marios Hadjieleftheriou, mhadji@gmail.com
 ******************************************************************************
 * Copyright (c) 2004, Marios Hadjieleftheriou
 *
 * All rights reserved.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
******************************************************************************/

#include <time.h>
#include <limits>
#include <spatialindex/tools/Tools.h>

#ifndef HAVE_SRAND48
#include <spatialindex/tools/rand48.h>
#endif

#include <cstring>

#if HAVE_PTHREAD_H
#if !defined(_POSIX_THREADS)
#include <sched.h>	// for sched_yield()
#endif
#endif

LibSupermesh_Tools::IndexOutOfBoundsException::IndexOutOfBoundsException(size_t i)
{
	std::ostringstream s;
	s << "Invalid index " << i;
	m_error = s.str();
}

std::string LibSupermesh_Tools::IndexOutOfBoundsException::what()
{
	return "IndexOutOfBoundsException: " + m_error;
}

LibSupermesh_Tools::IllegalArgumentException::IllegalArgumentException(std::string s) : m_error(s)
{
}

std::string LibSupermesh_Tools::IllegalArgumentException::what()
{
	return "IllegalArgumentException: " + m_error;
}

LibSupermesh_Tools::IllegalStateException::IllegalStateException(std::string s) : m_error(s)
{
}

std::string LibSupermesh_Tools::IllegalStateException::what()
{
	return "IllegalStateException: " + m_error;
}

LibSupermesh_Tools::EndOfStreamException::EndOfStreamException(std::string s) : m_error(s)
{
}

std::string LibSupermesh_Tools::EndOfStreamException::what()
{
	return "EndOfStreamException: " + m_error;
}

LibSupermesh_Tools::ResourceLockedException::ResourceLockedException(std::string s) : m_error(s)
{
}

std::string LibSupermesh_Tools::ResourceLockedException::what()
{
	return "ResourceLockedException: " + m_error;
}

LibSupermesh_Tools::NotSupportedException::NotSupportedException(std::string s) : m_error(s)
{
}

std::string LibSupermesh_Tools::NotSupportedException::what()
{
	return "NotSupportedException: " + m_error;
}

LibSupermesh_Tools::Variant::Variant() : m_varType(VT_EMPTY)
{
}

LibSupermesh_Tools::PropertySet::PropertySet(const byte* data)
{
	loadFromByteArray(data);
}

LibSupermesh_Tools::PropertySet::~PropertySet()
{
}

LibSupermesh_Tools::PropertySet::PropertySet()
{
}
void LibSupermesh_Tools::PropertySet::loadFromByteArray(const byte* ptr)
{
	m_propertySet.clear();

	uint32_t numberOfProperties;
	memcpy(&numberOfProperties, ptr, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	Variant v;

	for (uint32_t cIndex = 0; cIndex < numberOfProperties; ++cIndex)
	{
		std::string s(reinterpret_cast<const char*>(ptr));
		ptr += s.size() + 1;
		memcpy(&(v.m_varType), ptr, sizeof(VariantType));
		ptr += sizeof(VariantType);

		switch (v.m_varType)
		{
		case VT_SHORT:
			int16_t s;
			memcpy(&s, ptr, sizeof(int16_t));
			ptr += sizeof(int16_t);
			v.m_val.iVal = s;
			break;
		case VT_LONG:
			int32_t l;
			memcpy(&l, ptr, sizeof(int32_t));
			ptr += sizeof(int32_t);
			v.m_val.lVal = l;
			break;
		case VT_LONGLONG:
			int64_t ll;
			memcpy(&ll, ptr, sizeof(int64_t));
			ptr += sizeof(int64_t);
			v.m_val.llVal = ll;
			break;
		case VT_BYTE:
			byte b;
			memcpy(&b, ptr, sizeof(byte));
			ptr += sizeof(byte);
			v.m_val.bVal = b;
			break;
		case VT_FLOAT:
			float f;
			memcpy(&f, ptr, sizeof(float));
			ptr += sizeof(float);
			v.m_val.fltVal = f;
			break;
		case VT_DOUBLE:
			double d;
			memcpy(&d, ptr, sizeof(double));
			ptr += sizeof(double);
			v.m_val.dblVal = d;
			break;
		case VT_CHAR:
			char c;
			memcpy(&c, ptr, sizeof(char));
			ptr += sizeof(char);
			v.m_val.cVal = c;
			break;
		case VT_USHORT:
			uint16_t us;
			memcpy(&us, ptr, sizeof(uint16_t));
			ptr += sizeof(uint16_t);
			v.m_val.uiVal = us;
			break;
		case VT_ULONG:
			uint32_t ul;
			memcpy(&ul, ptr, sizeof(uint32_t));
			ptr += sizeof(uint32_t);
			v.m_val.ulVal = ul;
			break;
		case VT_ULONGLONG:
			uint64_t ull;
			memcpy(&ull, ptr, sizeof(uint64_t));
			ptr += sizeof(uint64_t);
			v.m_val.ullVal = ull;
			break;
		case VT_BOOL:
			byte bl;
			memcpy(&bl, ptr, sizeof(byte));
			ptr += sizeof(byte);
			v.m_val.blVal = (bl != 0);
			break;
		default:
			throw IllegalStateException(
				"Tools::PropertySet::PropertySet: Deserialization problem."
			);
		}

		m_propertySet.insert(std::pair<std::string, Variant>(s, v));
	}
}

uint32_t LibSupermesh_Tools::PropertySet::getByteArraySize()
{
	uint32_t size = sizeof(uint32_t);
	std::map<std::string, Variant>::iterator it;

	for (it = m_propertySet.begin(); it != m_propertySet.end(); ++it)
	{
		switch ((*it).second.m_varType)
		{
		case VT_SHORT:
			size += sizeof(int16_t);
			break;
		case VT_LONG:
			size += sizeof(int32_t);
			break;
		case VT_LONGLONG:
			size += sizeof(int64_t);
			break;
		case VT_BYTE:
			size += sizeof(byte);
			break;
		case VT_FLOAT:
			size += sizeof(float);
			break;
		case VT_DOUBLE:
			size += sizeof(double);
			break;
		case VT_CHAR:
			size += sizeof(char);
			break;
		case VT_USHORT:
			size += sizeof(uint16_t);
			break;
		case VT_ULONG:
			size += sizeof(uint32_t);
			break;
		case VT_ULONGLONG:
			size += sizeof(uint64_t);
			break;
		case VT_BOOL:
			size += sizeof(byte);
			break;
		default:
			throw NotSupportedException(
				"Tools::PropertySet::getSize: Unknown type."
			);
		}
		size += static_cast<uint32_t>((*it).first.size()) + 1 + sizeof(VariantType);
	}

	return size;
}

void LibSupermesh_Tools::PropertySet::storeToByteArray(byte** data, uint32_t& length)
{
	length = getByteArraySize();
	*data = new byte[length];
	byte* ptr = *data;

	uint32_t numberOfProperties = static_cast<uint32_t>(m_propertySet.size());
	memcpy(ptr, &numberOfProperties, sizeof(uint32_t));
	ptr += sizeof(uint32_t);

	std::map<std::string, Variant>::iterator it;

	for (it = m_propertySet.begin(); it != m_propertySet.end(); ++it)
	{
		size_t strSize = (*it).first.size();
		memcpy(ptr, (*it).first.c_str(), strSize);
		ptr += strSize;
		*ptr = 0;
		++ptr;

		memcpy(ptr, &((*it).second.m_varType), sizeof(VariantType));
		ptr += sizeof(VariantType);

		switch ((*it).second.m_varType)
		{
		case VT_SHORT:
			memcpy(ptr, &((*it).second.m_val.iVal), sizeof(int16_t));
			ptr += sizeof(int16_t);
			break;
		case VT_LONG:
			memcpy(ptr, &((*it).second.m_val.lVal), sizeof(int32_t));
			ptr += sizeof(int32_t);
			break;
		case VT_LONGLONG:
			memcpy(ptr, &((*it).second.m_val.llVal), sizeof(int64_t));
			ptr += sizeof(int64_t);
			break;
		case VT_BYTE:
			memcpy(ptr, &((*it).second.m_val.bVal), sizeof(byte));
			ptr += sizeof(byte);
			break;
		case VT_FLOAT:
			memcpy(ptr, &((*it).second.m_val.fltVal), sizeof(float));
			ptr += sizeof(float);
			break;
		case VT_DOUBLE:
			memcpy(ptr, &((*it).second.m_val.dblVal), sizeof(double));
			ptr += sizeof(double);
			break;
		case VT_CHAR:
			memcpy(ptr, &((*it).second.m_val.cVal), sizeof(char));
			ptr += sizeof(char);
			break;
		case VT_USHORT:
			memcpy(ptr, &((*it).second.m_val.uiVal), sizeof(uint16_t));
			ptr += sizeof(uint16_t);
			break;
		case VT_ULONG:
			memcpy(ptr, &((*it).second.m_val.ulVal), sizeof(uint32_t));
			ptr += sizeof(uint32_t);
			break;
		case VT_ULONGLONG:
			memcpy(ptr, &((*it).second.m_val.ullVal), sizeof(uint64_t));
			ptr += sizeof(uint64_t);
			break;
		case VT_BOOL:
			byte bl;
			bl = (*it).second.m_val.blVal;
			memcpy(ptr, &bl, sizeof(byte));
			ptr += sizeof(byte);
			break;
		default:
			throw NotSupportedException(
				"Tools::PropertySet::getData: Cannot serialize a variant of this type."
			);
		}
	}

	assert(ptr == (*data) + length);
}

LibSupermesh_Tools::Variant LibSupermesh_Tools::PropertySet::getProperty(std::string property)
{
   	std::map<std::string, Variant>::iterator it = m_propertySet.find(property);

   	if (it != m_propertySet.end()) return (*it).second;
   	else return Variant();
}

void LibSupermesh_Tools::PropertySet::setProperty(std::string property, Variant& v)
{
	std::pair<std::map<std::string, Variant>::iterator, bool> ret;
	std::map<std::string, Variant>::iterator it;

	ret = m_propertySet.insert(std::pair<std::string, Variant>(property, v));

	// If we weren't able to insert because it is already in the map
	// update our existing value
	if (ret.second == false) ret.first->second = v;
}

void LibSupermesh_Tools::PropertySet::removeProperty(std::string property)
{
   	std::map<std::string, Variant>::iterator it = m_propertySet.find(property);
   	if (it != m_propertySet.end()) m_propertySet.erase(it);
}

LibSupermesh_Tools::Interval::Interval() : m_type(IT_RIGHTOPEN), m_low(0.0), m_high(0.0)
{
}

LibSupermesh_Tools::Interval::Interval(IntervalType t, double l, double h) : m_type(t), m_low(l), m_high(h)
{
	assert(l < h);
}

LibSupermesh_Tools::Interval::Interval(double l, double h) : m_type(IT_RIGHTOPEN), m_low(l), m_high(h)
{
	assert(l < h);
}

LibSupermesh_Tools::Interval::Interval(const Interval& iv)
{
	m_low = iv.m_low;
	m_high = iv.m_high;
	m_type = iv.m_type;
}

LibSupermesh_Tools::IInterval& LibSupermesh_Tools::Interval::operator=(const LibSupermesh_Tools::IInterval& iv)
{
	if (this != &iv)
	{
		m_low = iv.getLowerBound();
		m_high = iv.getUpperBound();
		m_type = iv.getIntervalType();
	}

	return *this;
}

bool LibSupermesh_Tools::Interval::operator==(const Interval& iv) const
{
	if (
		m_type == iv.m_type &&
		m_low >= iv.m_low - std::numeric_limits<double>::epsilon() &&
		m_low <= iv.m_low + std::numeric_limits<double>::epsilon() &&
		m_high >= iv.m_high - std::numeric_limits<double>::epsilon() &&
		m_high <= iv.m_high + std::numeric_limits<double>::epsilon())
		return true;

	return false;
}

bool LibSupermesh_Tools::Interval::operator!=(const Interval& iv) const
{
	return ! (*this == iv);
}

double LibSupermesh_Tools::Interval::getLowerBound() const
{
	return m_low;
}

double LibSupermesh_Tools::Interval::getUpperBound() const
{
	return m_high;
}

void LibSupermesh_Tools::Interval::setBounds(double l, double h)
{
	assert(l <= h);

	m_low = l;
	m_high = h;
}

bool LibSupermesh_Tools::Interval::intersectsInterval(const IInterval& i) const
{
	return intersectsInterval(i.getIntervalType(), i.getLowerBound(), i.getUpperBound());
}

bool LibSupermesh_Tools::Interval::intersectsInterval(IntervalType type, const double low, const double high) const
{
	if (m_high < m_low)
		throw IllegalStateException(
			"Tools::Interval::intersectsInterval: high boundary is smaller than low boundary."
		);

	if (m_low > high || m_high < low) return false;
	if ((m_low > low && m_low < high) || (m_high > low && m_high < high)) return true;

	switch (m_type)
	{
	case IT_CLOSED:
		if (m_low == high)
		{
			if (type == IT_CLOSED || type == IT_LEFTOPEN) return true;
			else return false;
		}
		else if (m_high == low)
		{
			if (type == IT_CLOSED || type == IT_RIGHTOPEN) return true;
			else return false;
		}
		break;
	case IT_OPEN:
		if (m_low == high || m_high == low) return false;
		break;
	case IT_RIGHTOPEN:
		if (m_low == high)
		{
			if (type == IT_CLOSED || type == IT_LEFTOPEN) return true;
			else return false;
		}
		else if (m_high == low)
		{
			return false;
		}
		break;
	case IT_LEFTOPEN:
		if (m_low == high)
		{
			return false;
		}
		else if (m_high == low)
		{
			if (type == IT_CLOSED || type == IT_RIGHTOPEN) return true;
			else return false;
		}
		break;
	}

	return true;
}

bool LibSupermesh_Tools::Interval::containsInterval(const IInterval& i) const
{
	if (m_high < m_low)
		throw IllegalStateException(
			"Tools::Interval::containsInterval: high boundary is smaller than low boundary."
		);

	double low = i.getLowerBound();
	double high = i.getUpperBound();
	IntervalType type = i.getIntervalType();

	if (m_low < low && m_high > high) return true;
	if (m_low > low || m_high < high) return false;

	switch (m_type)
	{
	case IT_CLOSED:
		break;
	case IT_OPEN:
		if ((m_low == low && m_high == high && type != IT_OPEN) ||
			(m_low == low && (type == IT_CLOSED || type == IT_RIGHTOPEN)) ||
			(m_high == high && ( type == IT_CLOSED || type == IT_LEFTOPEN)))
			return false;
		break;
	case IT_RIGHTOPEN:
		if (m_high == high && (type == IT_CLOSED || type == IT_LEFTOPEN))
			return false;
		break;
	case IT_LEFTOPEN:
		if (m_low == low && (type == IT_CLOSED || type == IT_RIGHTOPEN))
			return false;
		break;
	}

	return true;
}

LibSupermesh_Tools::IntervalType LibSupermesh_Tools::Interval::getIntervalType() const
{
	return m_type;
}

LibSupermesh_Tools::Random::Random()
{
	m_pBuffer = 0;
	initDrand(static_cast<uint32_t>(time(0)), 0xD31A);
}

LibSupermesh_Tools::Random::Random(uint32_t seed, uint16_t xsubi0)
{
	m_pBuffer = 0;
	initDrand(seed, xsubi0);
}

LibSupermesh_Tools::Random::~Random()
{
	delete[] m_pBuffer;
}

void LibSupermesh_Tools::Random::initDrand(uint32_t seed, uint16_t xsubi0)
{
	m_pBuffer = new uint16_t[3];
	m_pBuffer[0] = static_cast<uint16_t>(xsubi0);
	uint32_t mask = 0xFFFF;
	m_pBuffer[1] = static_cast<uint16_t>(seed & mask);
	mask = mask << 16;
	m_pBuffer[2] = static_cast<uint16_t>((seed & mask) >> 16);

#ifdef BUILD_OS_CYGWIN
	srand48(*(reinterpret_cast<uint32_t*>(m_pBuffer)));
		// BUG: There is a bug in cygwin gcc 3.4.4. srand48 needs to be called
		// even if we are using the functions that take the seed as a parameter.
		// This does not affect random number generation, which still happens
		// using the seed provided as a parameter and not the one provided to srand48 :-S
#endif
}

int32_t LibSupermesh_Tools::Random::nextUniformLong()
{
	return jrand48(m_pBuffer);
}

uint32_t LibSupermesh_Tools::Random::nextUniformUnsignedLong()
{
	return static_cast<uint32_t>(nextUniformLong());
}

int32_t LibSupermesh_Tools::Random::nextUniformLong(int32_t low, int32_t high)
{
	return low + static_cast<int32_t>((high - low) * nextUniformDouble());
}

uint32_t LibSupermesh_Tools::Random::nextUniformUnsignedLong(uint32_t low, uint32_t high)
{
	return low + static_cast<uint32_t>((high - low) * nextUniformDouble());
}

int64_t LibSupermesh_Tools::Random::nextUniformLongLong()
{
	return static_cast<int64_t>(nextUniformUnsignedLongLong());
}

uint64_t LibSupermesh_Tools::Random::nextUniformUnsignedLongLong()
{
	uint64_t lh = static_cast<uint64_t>(nextUniformLong());
	uint64_t ll = static_cast<uint64_t>(nextUniformLong());
	uint64_t ret = (lh << 32) | ll;
	return ret;
}

int64_t LibSupermesh_Tools::Random::nextUniformLongLong(int64_t low, int64_t high)
{
	return low + static_cast<int64_t>((high - low) * nextUniformDouble());
}

uint64_t LibSupermesh_Tools::Random::nextUniformUnsignedLongLong(uint64_t low, uint64_t high)
{
	return low + static_cast<uint64_t>((high - low) * nextUniformDouble());
}

int16_t LibSupermesh_Tools::Random::nextUniformShort()
{
	return static_cast<int16_t>(nextUniformUnsignedShort());
}

uint16_t LibSupermesh_Tools::Random::nextUniformUnsignedShort()
{
	return nextUniformUnsignedLong() >> 16;
		// retain the high order bits.
}

double LibSupermesh_Tools::Random::nextUniformDouble()
{
	uint16_t* xsubi = reinterpret_cast<uint16_t*>(m_pBuffer);
	return erand48(xsubi);
}

double LibSupermesh_Tools::Random::nextUniformDouble(double low, double high)
{
	return (high - low) * nextUniformDouble() + low;
}

bool LibSupermesh_Tools::Random::flipCoin()
{
	if (nextUniformDouble() < 0.5) return true;
	return false;
}

#if HAVE_PTHREAD_H
LibSupermesh_Tools::LockGuard::LockGuard(pthread_mutex_t* pLock)
 : m_pLock(pLock)
{
	pthread_mutex_lock(m_pLock);
}

LibSupermesh_Tools::LockGuard::~LockGuard()
{
	pthread_mutex_unlock(m_pLock);
}
#endif

std::ostream& LibSupermesh_Tools::operator<<(std::ostream& os, const LibSupermesh_Tools::PropertySet& p)
{
	std::map<std::string, Variant>::const_iterator it;

	for (it = p.m_propertySet.begin(); it != p.m_propertySet.end(); ++it)
	{
		if (it != p.m_propertySet.begin()) os << ", ";

		switch ((*it).second.m_varType)
		{
		case VT_LONG:
			os << (*it).first << ": " << (*it).second.m_val.lVal;
			break;
		case VT_LONGLONG:
			os << (*it).first << ": " << (*it).second.m_val.llVal;
			break;
		case VT_BYTE:
			os << (*it).first << ": " << (*it).second.m_val.bVal;
			break;
		case VT_SHORT:
			os << (*it).first << ": " << (*it).second.m_val.iVal;
			break;
		case VT_FLOAT:
			os << (*it).first << ": " << (*it).second.m_val.fltVal;
			break;
		case VT_DOUBLE:
			os << (*it).first << ": " << (*it).second.m_val.dblVal;
			break;
		case VT_CHAR:
			os << (*it).first << ": " << (*it).second.m_val.cVal;
			break;
		case VT_USHORT:
			os << (*it).first << ": " << (*it).second.m_val.uiVal;
			break;
		case VT_ULONG:
			os << (*it).first << ": " << (*it).second.m_val.ulVal;
			break;
		case VT_ULONGLONG:
			os << (*it).first << ": " << (*it).second.m_val.ullVal;
			break;
		case VT_BOOL:
			os << (*it).first << ": " << (*it).second.m_val.blVal;
			break;
		case VT_PCHAR:
			os << (*it).first << ": " << (*it).second.m_val.pcVal;
			break;
		case VT_PVOID:
			os << (*it).first << ": ?";
			break;
		case VT_EMPTY:
			os << (*it).first << ": empty";
			break;
		default:
			os << (*it).first << ": unknown";
		}
	}

	return os;
}

std::ostream& LibSupermesh_Tools::operator<<(std::ostream& os, const LibSupermesh_Tools::Interval& iv)
{
	os << iv.m_type << " " << iv.m_low << " " << iv.m_high;
	return os;
}

//
// BufferedFile
//
LibSupermesh_Tools::BufferedFile::BufferedFile(uint32_t u32BufferSize)
: m_buffer(new char[u32BufferSize]), m_u32BufferSize(u32BufferSize), m_bEOF(true)
{
}

LibSupermesh_Tools::BufferedFile::~BufferedFile()
{
	m_file.close();
	delete[] m_buffer;
}

void LibSupermesh_Tools::BufferedFile::close()
{
	m_file.close();
}

bool LibSupermesh_Tools::BufferedFile::eof()
{
	return m_bEOF;
}

//
// BufferedFileReader
//
LibSupermesh_Tools::BufferedFileReader::BufferedFileReader()
{
}

LibSupermesh_Tools::BufferedFileReader::BufferedFileReader(const std::string& sFileName, uint32_t u32BufferSize)
: BufferedFile(u32BufferSize)
{
	open(sFileName);
}

void LibSupermesh_Tools::BufferedFileReader::open(const std::string& sFileName)
{
	m_bEOF = false;
	m_file.close(); m_file.clear();


	m_file.open(sFileName.c_str(), std::ios_base::in | std::ios_base::binary);
	if (! m_file.good())
		throw std::ios_base::failure("Tools::BufferedFileReader::BufferedFileReader: Cannot open file.");

	m_file.rdbuf()->pubsetbuf(m_buffer, m_u32BufferSize);
}

LibSupermesh_Tools::BufferedFileReader::~BufferedFileReader()
{
}

void LibSupermesh_Tools::BufferedFileReader::rewind()
{
	m_file.clear();
	m_file.seekg(0, std::ios_base::beg);
	if (! m_file.good())
		throw std::ios_base::failure("Tools::BufferedFileReader::rewind: seek failed.");

	m_bEOF = false;
}

void LibSupermesh_Tools::BufferedFileReader::seek(std::fstream::off_type offset)
{
	m_bEOF = false;
	m_file.clear();
	m_file.seekg(offset, std::ios_base::beg);
	if (! m_file.good())
		throw std::ios_base::failure("Tools::BufferedFileReader::seek: seek failed.");
}

uint8_t LibSupermesh_Tools::BufferedFileReader::readUInt8()
{
	if (m_bEOF) throw LibSupermesh_Tools::EndOfStreamException("");

	uint8_t ret;
	m_file.read(reinterpret_cast<char*>(&ret), sizeof(uint8_t));
	if (! m_file.good())
	{
		m_bEOF = true;
		throw LibSupermesh_Tools::EndOfStreamException("");
	}
	return ret;
}

uint16_t LibSupermesh_Tools::BufferedFileReader::readUInt16()
{
	if (m_bEOF) throw LibSupermesh_Tools::EndOfStreamException("");

	uint16_t ret;
	m_file.read(reinterpret_cast<char*>(&ret), sizeof(uint16_t));
	if (! m_file.good())
	{
		m_bEOF = true;
		throw LibSupermesh_Tools::EndOfStreamException("");
	}
	return ret;
}

uint32_t LibSupermesh_Tools::BufferedFileReader::readUInt32()
{
	if (m_bEOF) throw LibSupermesh_Tools::EndOfStreamException("");

	uint32_t ret;
	m_file.read(reinterpret_cast<char*>(&ret), sizeof(uint32_t));
	if (! m_file.good())
	{
		m_bEOF = true;
		throw LibSupermesh_Tools::EndOfStreamException("");
	}
	return ret;
}

uint64_t LibSupermesh_Tools::BufferedFileReader::readUInt64()
{
	if (m_bEOF) throw LibSupermesh_Tools::EndOfStreamException("");

	uint64_t ret;
	m_file.read(reinterpret_cast<char*>(&ret), sizeof(uint64_t));
	if (! m_file.good())
	{
		m_bEOF = true;
		throw LibSupermesh_Tools::EndOfStreamException("");
	}
	return ret;
}

float LibSupermesh_Tools::BufferedFileReader::readFloat()
{
	if (m_bEOF) throw LibSupermesh_Tools::EndOfStreamException("");

	float ret;
	m_file.read(reinterpret_cast<char*>(&ret), sizeof(float));
	if (! m_file.good())
	{
		m_bEOF = true;
		throw LibSupermesh_Tools::EndOfStreamException("");
	}
	return ret;
}

double LibSupermesh_Tools::BufferedFileReader::readDouble()
{
	if (m_bEOF) throw LibSupermesh_Tools::EndOfStreamException("");

	double ret;
	m_file.read(reinterpret_cast<char*>(&ret), sizeof(double));
	if (! m_file.good())
	{
		m_bEOF = true;
		throw LibSupermesh_Tools::EndOfStreamException("");
	}
	return ret;
}

bool LibSupermesh_Tools::BufferedFileReader::readBoolean()
{
	if (m_bEOF) throw LibSupermesh_Tools::EndOfStreamException("");

	bool ret;
	m_file.read(reinterpret_cast<char*>(&ret), sizeof(bool));
	if (! m_file.good())
	{
		m_bEOF = true;
		throw LibSupermesh_Tools::EndOfStreamException("");
	}
	return ret;
}

std::string LibSupermesh_Tools::BufferedFileReader::readString()
{
	if (m_bEOF) throw LibSupermesh_Tools::EndOfStreamException("");

	uint32_t len;
	m_file.read(reinterpret_cast<char*>(&len), sizeof(uint32_t));
	if (! m_file.good())
	{
		m_bEOF = true;
		throw LibSupermesh_Tools::EndOfStreamException("");
	}

	std::string::value_type* buf = new std::string::value_type[len];
	m_file.read(reinterpret_cast<char*>(buf), len * sizeof(std::string::value_type));
	if (! m_file.good())
	{
		delete[] buf;
		m_bEOF = true;
		throw LibSupermesh_Tools::EndOfStreamException("");
	}

	std::string ret(buf, len);
	delete[] buf;

	return ret;
}

void LibSupermesh_Tools::BufferedFileReader::readBytes(uint32_t u32Len, byte** pData)
{
	if (m_bEOF) throw LibSupermesh_Tools::EndOfStreamException("");

	*pData = new byte[u32Len];
	m_file.read(reinterpret_cast<char*>(*pData), u32Len);
	if (! m_file.good())
	{
		delete[] *pData;
		m_bEOF = true;
		throw LibSupermesh_Tools::EndOfStreamException("");
	}
}

//
// BufferedFileWriter
//
LibSupermesh_Tools::BufferedFileWriter::BufferedFileWriter()
{
	open("");
}

LibSupermesh_Tools::BufferedFileWriter::BufferedFileWriter(const std::string& sFileName, FileMode mode, uint32_t u32BufferSize)
: BufferedFile(u32BufferSize)
{
	open(sFileName, mode);
}

LibSupermesh_Tools::BufferedFileWriter::~BufferedFileWriter()
{
	m_file.flush();
}

void LibSupermesh_Tools::BufferedFileWriter::open(const std::string& sFileName, FileMode mode)
{
	m_bEOF = false;
	m_file.close(); m_file.clear();

	if (mode == CREATE)
	{
		m_file.open(sFileName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
		if (! m_file.good())
			throw std::ios_base::failure("Tools::BufferedFileWriter::open: Cannot open file.");
	}
	else if (mode == APPEND)
	{
		// Idiotic fstream::open truncates an existing file anyway, if it is only opened
		// for output (no ios_base::in flag)!! On the other hand, if a file does not exist
		// and the ios_base::in flag is specified, then the open fails!!

		m_file.open(sFileName.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
		if (! m_file.good())
		{
			m_file.clear();
			m_file.open(sFileName.c_str(), std::ios_base::out | std::ios_base::binary);
			if (! m_file.good())
				throw std::ios_base::failure("Tools::BufferedFileWriter::open: Cannot open file.");
		}
		else
		{
			m_file.seekp(0, std::ios_base::end);
			if (! m_file.good())
				throw std::ios_base::failure("Tools::BufferedFileWriter::open: Cannot open file.");
		}
	}
	else
		throw LibSupermesh_Tools::IllegalArgumentException("Tools::BufferedFileWriter::open: Unknown mode.");
}

void LibSupermesh_Tools::BufferedFileWriter::rewind()
{
	m_bEOF = false;
	m_file.clear();
	m_file.seekp(0, std::ios_base::beg);
	if (! m_file.good())
		throw std::ios_base::failure("Tools::BufferedFileWriter::rewind: seek failed.");
}

void LibSupermesh_Tools::BufferedFileWriter::seek(std::fstream::off_type offset)
{
	m_bEOF = false;
	m_file.clear();
	m_file.seekp(offset, std::ios_base::beg);
	if (! m_file.good())
		throw std::ios_base::failure("Tools::BufferedFileWriter::seek: seek failed.");
}

void LibSupermesh_Tools::BufferedFileWriter::write(uint8_t i)
{
	m_file.write(reinterpret_cast<const char*>(&i), sizeof(uint8_t));
	if (! m_file.good()) throw std::ios_base::failure("");
}

void LibSupermesh_Tools::BufferedFileWriter::write(uint16_t i)
{
	m_file.write(reinterpret_cast<const char*>(&i), sizeof(uint16_t));
	if (! m_file.good()) throw std::ios_base::failure("");
}

void LibSupermesh_Tools::BufferedFileWriter::write(uint32_t i)
{
	m_file.write(reinterpret_cast<const char*>(&i), sizeof(uint32_t));
	if (! m_file.good()) throw std::ios_base::failure("");
}

void LibSupermesh_Tools::BufferedFileWriter::write(uint64_t i)
{
	m_file.write(reinterpret_cast<const char*>(&i), sizeof(uint64_t));
	if (! m_file.good()) throw std::ios_base::failure("");
}

void LibSupermesh_Tools::BufferedFileWriter::write(float i)
{
	m_file.write(reinterpret_cast<const char*>(&i), sizeof(float));
	if (! m_file.good()) throw std::ios_base::failure("");
}

void LibSupermesh_Tools::BufferedFileWriter::write(double i)
{
	m_file.write(reinterpret_cast<const char*>(&i), sizeof(double));
	if (! m_file.good()) throw std::ios_base::failure("");
}

void LibSupermesh_Tools::BufferedFileWriter::write(bool b)
{
	m_file.write(reinterpret_cast<const char*>(&b), sizeof(bool));
	if (! m_file.good()) throw std::ios_base::failure("");
}

void LibSupermesh_Tools::BufferedFileWriter::write(const std::string& s)
{
	uint32_t len = static_cast<uint32_t>(s.size());
	m_file.write(reinterpret_cast<const char*>(&len), sizeof(uint32_t));
	if (! m_file.good()) throw std::ios_base::failure("");
	m_file.write(reinterpret_cast<const char*>(s.c_str()), len * sizeof(std::string::value_type));
	if (! m_file.good()) throw std::ios_base::failure("");
}

void LibSupermesh_Tools::BufferedFileWriter::write(uint32_t u32Len, byte* pData)
{
	m_file.write(reinterpret_cast<const char*>(pData), u32Len);
	if (! m_file.good()) throw std::ios_base::failure("");
}

//
// TemporaryFile
//
LibSupermesh_Tools::TemporaryFile::TemporaryFile()
{

#ifdef _MSC_VER

#ifndef L_tmpnam_s
// MSVC 2003 doesn't have tmpnam_s, so we'll have to use the old functions

    char* tmpName = NULL;
    tmpName = tmpnam( NULL );

    if (tmpName == NULL)
        throw std::ios_base::failure("Tools::TemporaryFile: Cannot create temporary file name.");

#else
	char tmpName[L_tmpnam_s];
	errno_t err = tmpnam_s(tmpName, L_tmpnam_s);
	if (err)
		throw std::ios_base::failure("Tools::TemporaryFile: Cannot create temporary file name.");

#endif
	if (tmpName[0] == '\\')
		m_sFile = std::string(tmpName + 1);
	else
		m_sFile = std::string(tmpName);

#else
	char tmpName[7] = "XXXXXX";
	if (mktemp(tmpName) == 0)
		throw std::ios_base::failure("Tools::TemporaryFile: Cannot create temporary file name.");
	m_sFile = tmpName;
#endif

	m_pFile = new LibSupermesh_Tools::BufferedFileWriter(m_sFile, LibSupermesh_Tools::CREATE);
}

LibSupermesh_Tools::TemporaryFile::~TemporaryFile()
{
	delete m_pFile;

#ifdef _MSC_VER
	_unlink(m_sFile.c_str());
#else
	std::remove(m_sFile.c_str());
#endif
}

void LibSupermesh_Tools::TemporaryFile::rewindForReading()
{
	LibSupermesh_Tools::BufferedFileReader* br = dynamic_cast<LibSupermesh_Tools::BufferedFileReader*>(m_pFile);
	if (br != 0)
		m_pFile->rewind();
	else
	{
		delete m_pFile;
		m_pFile = new LibSupermesh_Tools::BufferedFileReader(m_sFile);
	}
}

void LibSupermesh_Tools::TemporaryFile::rewindForWriting()
{
	LibSupermesh_Tools::BufferedFileWriter* bw = dynamic_cast<LibSupermesh_Tools::BufferedFileWriter*>(m_pFile);
	if (bw != 0)
		m_pFile->rewind();
	else
	{
		delete m_pFile;
		m_pFile = new LibSupermesh_Tools::BufferedFileWriter(m_sFile);
	}
}

bool LibSupermesh_Tools::TemporaryFile::eof()
{
	return m_pFile->eof();
}

std::string LibSupermesh_Tools::TemporaryFile::getFileName() const
{
	return m_sFile;
}

uint8_t LibSupermesh_Tools::TemporaryFile::readUInt8()
{
	LibSupermesh_Tools::BufferedFileReader* br = dynamic_cast<LibSupermesh_Tools::BufferedFileReader*>(m_pFile);
	if (br == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::readUInt8: file not open for reading.");

	return br->readUInt8();
}

uint16_t LibSupermesh_Tools::TemporaryFile::readUInt16()
{
	LibSupermesh_Tools::BufferedFileReader* br = dynamic_cast<LibSupermesh_Tools::BufferedFileReader*>(m_pFile);
	if (br == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::readUInt16: file not open for reading.");

	return br->readUInt16();
}

uint32_t LibSupermesh_Tools::TemporaryFile::readUInt32()
{
	LibSupermesh_Tools::BufferedFileReader* br = dynamic_cast<LibSupermesh_Tools::BufferedFileReader*>(m_pFile);
	if (br == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::readUInt32: file not open for reading.");

	return br->readUInt32();
}

uint64_t LibSupermesh_Tools::TemporaryFile::readUInt64()
{
	LibSupermesh_Tools::BufferedFileReader* br = dynamic_cast<LibSupermesh_Tools::BufferedFileReader*>(m_pFile);
	if (br == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::readUInt64: file not open for reading.");

	return br->readUInt64();
}

float LibSupermesh_Tools::TemporaryFile::readFloat()
{
	LibSupermesh_Tools::BufferedFileReader* br = dynamic_cast<LibSupermesh_Tools::BufferedFileReader*>(m_pFile);
	if (br == 0)
		throw std::ios_base::failure("LibSupermesh_Tools::TemporaryFile::readFloat: file not open for reading.");

	return br->readFloat();
}

double LibSupermesh_Tools::TemporaryFile::readDouble()
{
	LibSupermesh_Tools::BufferedFileReader* br = dynamic_cast<LibSupermesh_Tools::BufferedFileReader*>(m_pFile);
	if (br == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::readDouble: file not open for reading.");

	return br->readDouble();
}

std::string LibSupermesh_Tools::TemporaryFile::readString()
{
	LibSupermesh_Tools::BufferedFileReader* br = dynamic_cast<LibSupermesh_Tools::BufferedFileReader*>(m_pFile);
	if (br == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::readString: file not open for reading.");

	return br->readString();
}

void LibSupermesh_Tools::TemporaryFile::readBytes(uint32_t u32Len, byte** pData)
{
	LibSupermesh_Tools::BufferedFileReader* br = dynamic_cast<LibSupermesh_Tools::BufferedFileReader*>(m_pFile);
	if (br == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::readString: file not open for reading.");

	return br->readBytes(u32Len, pData);
}

void LibSupermesh_Tools::TemporaryFile::write(uint8_t i)
{
	LibSupermesh_Tools::BufferedFileWriter* bw = dynamic_cast<LibSupermesh_Tools::BufferedFileWriter*>(m_pFile);
	if (bw == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::write: file not open for writing.");

	return bw->write(i);
}

void LibSupermesh_Tools::TemporaryFile::write(uint16_t i)
{
	LibSupermesh_Tools::BufferedFileWriter* bw = dynamic_cast<LibSupermesh_Tools::BufferedFileWriter*>(m_pFile);
	if (bw == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::write: file not open for writing.");

	return bw->write(i);
}

void LibSupermesh_Tools::TemporaryFile::write(uint32_t i)
{
	LibSupermesh_Tools::BufferedFileWriter* bw = dynamic_cast<LibSupermesh_Tools::BufferedFileWriter*>(m_pFile);
	if (bw == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::write: file not open for writing.");

	return bw->write(i);
}

void LibSupermesh_Tools::TemporaryFile::write(uint64_t i)
{
	LibSupermesh_Tools::BufferedFileWriter* bw = dynamic_cast<LibSupermesh_Tools::BufferedFileWriter*>(m_pFile);
	if (bw == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::write: file not open for writing.");

	return bw->write(i);
}

void LibSupermesh_Tools::TemporaryFile::write(float i)
{
	LibSupermesh_Tools::BufferedFileWriter* bw = dynamic_cast<LibSupermesh_Tools::BufferedFileWriter*>(m_pFile);
	if (bw == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::write: file not open for writing.");

	return bw->write(i);
}

void LibSupermesh_Tools::TemporaryFile::write(double i)
{
	LibSupermesh_Tools::BufferedFileWriter* bw = dynamic_cast<LibSupermesh_Tools::BufferedFileWriter*>(m_pFile);
	if (bw == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::write: file not open for writing.");

	return bw->write(i);
}

void LibSupermesh_Tools::TemporaryFile::write(const std::string& s)
{
	LibSupermesh_Tools::BufferedFileWriter* bw = dynamic_cast<LibSupermesh_Tools::BufferedFileWriter*>(m_pFile);
	if (bw == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::write: file not open for writing.");

	return bw->write(s);
}

void LibSupermesh_Tools::TemporaryFile::write(uint32_t u32Len, byte* pData)
{
	LibSupermesh_Tools::BufferedFileWriter* bw = dynamic_cast<LibSupermesh_Tools::BufferedFileWriter*>(m_pFile);
	if (bw == 0)
		throw std::ios_base::failure("Tools::TemporaryFile::write: file not open for writing.");

	return bw->write(u32Len, pData);
}
