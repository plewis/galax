//
//  split.hpp
//  galax
//
//  Created by Paul O. Lewis on 10/25/14.
//  Copyright (c) 2014 Paul O. Lewis. All rights reserved.
//

#ifndef GALAX_SPLIT_HPP
#define GALAX_SPLIT_HPP

#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>
#include <vector>
#include <set>
#include <sstream>
#include "xgalax.hpp"

namespace galax
{

class Split
	{	
    public:
        typedef unsigned long           split_t;        //!< type for splits, Split::InvertSplit assumes this is an unsigned integer type
        typedef std::vector<split_t>    SplitTVect;     //!< vector of split_t
        typedef std::set<Split>         SplitSet;       //!< set of Split objects

                                        Split();
							            Split(const Split & other);
							            ~Split();

        void                            copy(const Split & other);
		
		std::string 	                getDimensionInfo();
		unsigned 		                getNTaxa() const;
		void                            setNTaxa(unsigned n);
		void 			                calcNUnits(unsigned ntax);

		unsigned			            calcComplexity() const;

        std::string                     writeToString() const;
        void                            readFromString(const std::string s);
        void                            readFromStream(std::istream & in);

        void                            createFromPattern(std::string s);
		void	 			            createAndAppendPatternRepresentation(std::string &) const;
		std::string			            createIdRepresentation() const;
		std::string			            createPatternRepresentation() const;
        std::string                     createNewickRepresentation(bool zero_based = false) const;

        unsigned			            countOnBits() const;
		unsigned			            countOffBits() const;
		int 				            cmp(const Split & other) const;
		bool 				            equals(const Split & other) const;
		bool 				            isBitSet(unsigned t) const;
		bool				            isCompatible(const Split & other) const;
		bool 				            isLessThan(const Split & other) const;
		bool 				            subsumedIn(const Split & other, unsigned startUnit = 0) const;

        std::vector<unsigned>           getOnList() const;
        std::vector<unsigned>           getOffList() const;
        void                            setExcluded(const std::vector<unsigned> excl);
        std::vector<unsigned>           getExcludedList() const;

        void                            setOnSymbol(char c);
        void                            setOffSymbol(char c);
        void                            setExcludedSymbol(char c);

        char                            getOnSymbol() const;
        char                            getOffSymbol() const;
        char                            getExcludedSymbol() const;

		bool				            operator!=(const Split & other) const;
		Split 				            operator&(const Split & other) const;
		Split &                         operator&=(const Split & other);
		Split 				            operator^(const Split & other) const;
		Split &                         operator^=(const Split & other);
		Split 				            operator|(const Split & other) const;
		Split &                         operator|=(const Split & other);
		bool				            operator<(const Split & other) const;
		bool				            operator==(const Split & other) const;
		Split &                         operator=(const Split & other);

		//temporary! friend std::istream &           operator>>(std::istream & in, Split & s);

        //@POL 13-Nov-2007 it is currently a mystery why uncommenting the line below leads to hundreds of VC error C2679 messages: 
        // binary '<<' : no operator found which takes a right-hand operand of type <whatever> (or there is no acceptable conversion)
        //friend std::string &    operator<<(std::string & out, const Split & s);

		void 				            clear();
        void                            reset();
        void                            resetNUnits(unsigned ntax);

		void 				            combineWith(const Split & other);
		void 				            intersectWith(const Split & other);

		void 				            setBit(unsigned t);
        void 				            setBits(const std::vector<unsigned> bits_to_set);
		void 				            unsetBit(unsigned t);
        void 				            unsetBits(const std::vector<unsigned> bits_to_unset);

		void 				            invertSplit();
			
	private:

		void 				            resize();
        void                            getOnListImpl(std::vector<unsigned> & v) const;
        void                            getOffListImpl(std::vector<unsigned> & v) const;
        void                            writeImpl(std::string & out) const;
		
		SplitTVect                      _unit;              //!< is the vector of split units
		const unsigned		            _bits_per_unit;     //!< is the number of bits in a variable of type split_t
		const split_t		            _split_unity;       //!< is a split_t variable with only the least significant bit set
		unsigned		                _split_ntax;		//!< is the number of taxa currently under consideration
		unsigned		                _nunits;			//!< is the length of the array necessary to represent a split
		char			                _on_symbol;         //!< is the symbol used to represent bits that are set (i.e. "on")
		char			                _off_symbol;		//!< is the symbol used to represent bits that have been cleared (i.e. "off")
		char			                _excl_symbol;       //!< is the symbol used to represent bits that have been excluded (i.e. should not be considered off or on)
        std::vector<unsigned>           _excl_bits;         //!< is a sorted vector containing bit positions corresponding to excluded taxa (lowest possible value is 0)
        
	public:
		typedef boost::shared_ptr< Split > SplitShPtr;
	};

/**
*   Stream input operator used to read in a Split object from an istream.
*/
//temporary! std::istream & operator>>(std::istream & in, Split & s);

/**
*   Fills \c missing with splits that are in \c ref_tree but absent in \c test_tree.
*/
unsigned findSplitsAbsentInTestTree(const Split::SplitSet & ref_tree, const Split::SplitSet & test_tree, Split::SplitSet & missing);

/**
*   Calls Clear().
*/
inline Split::Split()
  : _bits_per_unit((CHAR_BIT)*sizeof(split_t)), _split_unity((split_t)1)
	{
	clear();
	}

/**
*   Sets #_split_ntax, #_nunits, #_on_symbol, #_off_symbol, #_excl_symbol to default values and
*   resizes the #_units vector to have a length of #_nunits and fills it with 0s.
*/
inline void Split::clear()
	{
    _split_ntax         = 4;
    _nunits             = 1;
    _on_symbol          = '*';
    _off_symbol         = '-';
    _excl_symbol		= 'x';
    _excl_bits.clear();
	_unit.resize(_nunits, (split_t)0);
    std::fill(_unit.begin(), _unit.end(), (split_t)0);
	}

/**
*   Sets each bit in the #_units vector to 0 (unset state), but does not change anything else. Useful for cleaning the
*   slate without affecting the dimensions (that is, data members #_split_ntax and #_nunits are not changed).
*/
inline void Split::reset()
	{
    std::fill(_unit.begin(), _unit.end(), (split_t)0);
	}

/**
*   Calls calcNUnits(), then resize() and finally reset().
*/
inline void Split::resetNUnits(unsigned ntax)
	{
    calcNUnits(ntax);
    resize();
    reset();
	}

/**
*   Creates a copy of the supplied \c other Split object by calling operator=.
*/
inline Split::Split(const Split & other)
  : _bits_per_unit((CHAR_BIT)*sizeof(split_t)), _split_unity((split_t)1)
	{
	*this = other;
	}

/**
*   Creates a copy of the supplied \c other Split object.
*/
inline void Split::copy(const Split & other)
	{
    _split_ntax		= other._split_ntax;
    _nunits			= other._nunits;
    _on_symbol		= other._on_symbol;
    _off_symbol		= other._off_symbol;
    _excl_symbol	= other._excl_symbol;
	_excl_bits.resize(other._excl_bits.size(), (split_t)0);
	std::copy(other._excl_bits.begin(), other._excl_bits.end(), _excl_bits.begin());
	_unit.resize(_nunits, (split_t)0);
	std::copy(other._unit.begin(), other._unit.end(), _unit.begin());
	}

/**
*   Creates a copy of the Split object \c other by calling Split::copy, then returns reference to *this.
*/
inline Split & Split::operator=(const Split & other)
	{
    copy(other);
	return *this;
	}

/**
*   Destructor does nothing because all data members are automatically deleted.
*/
inline Split::~Split()
	{
	}

/**
*   Sets non-constant data members using a pattern supplied in the form of the string \c s consisting of a sequence of
*   #_on_symbol, #_off_symbol and (possibly) #_excl_symbol characters. The data members #_on_symbol, #_off_symbol and
*   #_excl_symbol are not modified. For example, assuming #_bits_per_unit = 8, calling the function with the pattern 
*   \c --*---***- sets #_split_ntax to 10, #_nunits to 2,  #_unit[0] to 142 (binary 10001110), and #_unit[1] to 0.
*/
inline void Split::createFromPattern(std::string s)
    {
    unsigned slen = (unsigned)s.size();
    std::vector<unsigned> on_bits;
	_excl_bits.clear();
    for (unsigned k = 0; k < slen; ++k)
        {
        if (s[k] == _on_symbol)
            on_bits.push_back(k);
        else if (s[k] == _excl_symbol)
            _excl_bits.push_back(k);
        else if (s[k] != _off_symbol)
            throw XGalax(str(boost::format("character in pattern (%c) not recognized as either the on symbol (%c), off symbol (%c) or excluded symbol (%c)") % s[k] % _on_symbol % _off_symbol % _excl_symbol));
        }

    // No funny characters were found in the supplied pattern string, so we have a green light to build the split
    calcNUnits(slen);   // sets split_ntax and nunits
	_unit.resize(_nunits, (split_t)0);
    for (std::vector<unsigned>::const_iterator it = on_bits.begin(); it != on_bits.end(); ++it)
        {
        unsigned k = *it;
        unsigned i = k/_bits_per_unit;
        unsigned j = k % _bits_per_unit;
        _unit[i] |= (_split_unity << j);
        }
    }

/**
*   Recomputes those data members (#_split_ntax and #_nunits) that depend on the number of taxa in the tree. The data
*   member #_split_ntax is set to the specified number of taxa \c ntax. The data member #_nunits equals the number of
*   units (each holding #_bits_per_unit bits) needed to accommodate #_split_ntax taxa.
*/
inline void Split::calcNUnits(unsigned ntax)
	{
	if (ntax == 0)
		{
	    _split_ntax = 0;
		_nunits = 1;
		}
    else
        {
	    _split_ntax = ntax;

        // Suppose bits_per_unit = 8
        //   7 taxa requires 1 unit  ==> 1 + (7-1)/8 = 1 + 6/8 = 1 + 0 = 1
        //   8 taxa requires 1 unit  ==> 1 + (8-1)/8 = 1 + 7/8 = 1 + 0 = 1
        //   9 taxa requires 2 units ==> 1 + (9-1)/8 = 1 + 8/8 = 1 + 1 = 2
	    _nunits = 1 + ((_split_ntax - 1)/_bits_per_unit);
       }
	}

/**
*   Serializes a Split object to the supplied string \c out. Split objects output using this method can be input using 
*   the corresponding function readFromStream(). The information is space-delimited, and in the following order: #_bits_per_unit,
*   #_split_ntax, #_on_symbol, #_off_symbol, #_excl_symbol, pattern representing #_unit vector. #_nunits is not included 
*   because it can be recalculated from knowledge of #_bits_per_unit and #_split_ntax using the function calcNUnits(). 
*   Likewise, the #_unit and #_excl_bits vectors are not output explicitly because they can be constructed from the 
*   pattern that is output (and the pattern has the addition benefit of being human-readable). Note that \c out is 
*   cleared before any values are appended.
*/
inline void Split::writeImpl(std::string & out) const
	{
    out = boost::str(boost::format("%c%c%c") % _on_symbol % _off_symbol % _excl_symbol);
    out += createPatternRepresentation();
	}

/**
*   Calls writeImpl(), then returns the resulting string.
*/
inline std::string Split::writeToString() const
    {
    std::string out;
    writeImpl(out);
	return out;
    }

/**
*   Inputs from supplied stream \c in all the information needed to fully construct a Split object. Use the
*   corresponding writeImpl(), writeToString(), or operator<<() functions to save splits in a format that can be
*   read in using this operator. On failure, an XGalax exception will be thrown and \c s will be cleared
*   (that is, returned to a state identical to a default-constructed Split object).
*/
inline void Split::readFromStream(std::istream & in)
	{
    std::string s;
    in >> s;

    // If the above read operation failed, bail out now
    if (!in)
        throw XGalax("problem reading split from input stream");

    if (s.size() < 3)
        throw XGalax(boost::str(boost::format("split definition read from input stream (%s) was too short (that is, at minimum, the on, off and excluded symbols must be present)") % s));
    _on_symbol = s[0];
    _off_symbol = s[1];
    _excl_symbol = s[2];

    s.erase(s.begin(), s.begin() + 3);

    createFromPattern(s);
	}

/**
*   Creates an input string stream from \c in, then calls readFromStream().
*/
inline void Split::readFromString(const std::string in)
    {
    std::istringstream inss(in);
    readFromStream(inss);
    }

/**
*   Inputs from supplied ifstream \c in all the information needed to fully construct the Split \c s. Use operator<<() to
*   output Split objects so that they can be read in using this operator. On failure, an XGalax exception will be
*   thrown and \c s will be cleared (returned to a state identical to a default-constructed Split object).
*/
//temporary! std::istream & operator>>(std::istream & in, Split & s)
//	{
//    s.readFromStream(in);
//    return in;
//	}

/**
*   Creates a coded representation in which #_on_symbol represents bits that are "on" (that is, set to 1), #_off_symbol
*   represents "off" bits (that is, set to 0), and #_excl_symbol represent bits that have been excluded. Uses function
*   createAndAppendPatternRepresentation() to do the actual work.
*/
inline std::string Split::createPatternRepresentation() const
	{
	std::string s;
	s.reserve(_split_ntax);
	createAndAppendPatternRepresentation(s);
	return s;
	}

/**
*   Creates a coded representation in which #_on_symbol represents bits that are "on" (that is, set to 1), #_off_symbol
*   represents "off" bits (that is, set to 0), and #_excl_symbol represent bits that have been excluded. This function is
*   employed by createPatternRepresentation(). Note that taxon 0 occupies the leftmost position in the resulting
*   string, whereas in underlying representation things are not so tidy. Here is an example involving 10 taxa where 
*   #_bits_per_unit = 8 and taxa 1, 7 and 9 are "on":
*
*       +--------+--------+
*       |10000010|00000010|
*       +--------+--------+
*        _unit[0] _unit[1]
*
*   This split would be output as follows, assuming #_on_symbol is an asterisk and #_off_symbol is a hyphen:
*
*               -*-----*-*
*       taxon 0 ^        ^ taxon 9
*
*   Note that "taxon 9" actually represents the 10th. taxon.
*/
inline void Split::createAndAppendPatternRepresentation(std::string & s) const
	{
	if (_split_ntax > 0)
		{
		unsigned ntax_added = 0;
		for (unsigned i = 0; i < _nunits; ++i)
			{
			for (unsigned j = 0; j < _bits_per_unit; ++j)
				{
                split_t bitmask = (_split_unity << j);
				bool bit_is_set = ((_unit[i] & bitmask) > (split_t)0);
                if (!_excl_bits.empty() && std::binary_search(_excl_bits.begin(), _excl_bits.end(), ntax_added))
                    s += _excl_symbol;
                else if (bit_is_set)
                    s += _on_symbol;
                else
                    s += _off_symbol;
				if (++ntax_added == _split_ntax)
					return;
				}
			}
		}
	}

/**
*   Creates a newick tree description representing the split and returns as a string. Each taxon is represented by an
*   integer starting at 0 (if #_zero_based is true) or 1 (if #_zero_based is false). The taxa representing "on" bits are
*   listed first in the tree description. Thus, if bits 1, 7 and 9 were "on", and 0, 2-6, and 8 were "off", the tree 
*   description would look like this if #_zero_based was false:
*
*       (2,8,10,(1,3,4,5,6,7,9))
*
*   Note that no terminating semicolon is added by this function.
*/
inline std::string Split::createNewickRepresentation(bool zero_based) const
	{
    unsigned taxon_number, taxon_index;
	std::string s;
	s += '(';

    // Spit out "on" bits first
    std::vector<unsigned> on_bits;
    getOnListImpl(on_bits);
    for (std::vector<unsigned>::const_iterator i = on_bits.begin(); i != on_bits.end(); ++i)
        {
        taxon_index = *i;
        taxon_number = (zero_based ? taxon_index : taxon_index + 1);
        s += str(boost::format("%d,") % taxon_number);
        }

    // Now for the "off" bits
    std::vector<unsigned> off_bits;
    getOffListImpl(off_bits);
	if (off_bits.size() > 0)
		{
		s += '(';
		std::vector<unsigned>::const_iterator j = off_bits.begin();
		taxon_index = *j;
		taxon_number = (zero_based ? taxon_index : taxon_index + 1);
		s += str(boost::format("%d") % taxon_number);
		for (++j; j != off_bits.end(); ++j)
			{
			taxon_index = *j;
			taxon_number = (zero_based ? taxon_index : taxon_index + 1);
			s += str(boost::format(",%d") % taxon_number);
			}
		s += ')';
		}
	s += ')';
	return s;
	}

/**
*   Creates a representation of the Split in which each unit is simply output as a decimal value. This representation is
*   not easily interpreted, as each unit comprises bits for several taxa (and it is the binary representation that makes
*   sense). The units are output with #_unit[0] on the left. Here is an example involving #_bits_per_unit = 8:
*
*       +--------+--------+
*       |10000010|00000010|
*       +--------+--------+
*        _unit[0] _unit[1]
*
*   The string produced from the above #_unit vector would be:
*
*       130 2
*
*/
inline std::string Split::createIdRepresentation() const
	{
	std::string s;
    assert(_unit.size() > 0);
    s += boost::str(boost::format("%d") % _unit[0]);
	for (unsigned i = 1; i < _nunits; ++i)
        s += boost::str(boost::format(" %d") % _unit[i]);
	return s;
	}

/**
*   Returns the number of bits that are set to one (that is, are "on") in the Split.
*/
inline unsigned Split::countOnBits() const
	{
    // This assumes that the unused bits in the last unit are all 0.
    unsigned num_bits_set = 0;

    for (SplitTVect::const_iterator i = _unit.begin(); i != _unit.end(); ++i)
        {
        // The  Kernighan method of counting bits is used below: see exercise 2-9 in the Kernighan and Ritchie book.
        // As an example, consider v = 10100010
        // c = 0:
        //   v     = 10100010
        //   v - 1 = 10100001
        //   ----------------
        //   new v = 10100000
        // c = 1:
        //   v     = 10100000
        //   v - 1 = 10011111
        //   ----------------
        //   new v = 10000000
        // c = 2:
        //   v     = 10000000
        //   v - 1 = 01111111
        //   ----------------
        //   new v = 00000000
        // c = 3:
        //   break out of loop because v = 0
        // 
        split_t v = *i;
        unsigned c = 0;
        for (; v; ++c)
            {
            v &= v - 1;
            }
        num_bits_set += c;
        }
    return num_bits_set;
    }
	
/**
*   Returns the number of bits that are set to zero (that is, are "off") in the Split. This function is slower than the
*   function countOnBits unless #_excl_bits is empty, in which case the two functions will be essentially equal in speed.
*   The problem is that countOnBits does not need to check for excluded bits because excluded bits are guaranteed to be
*   off. If the #_excl_bits vector is empty, then countOffBits is slower only because it needs to invert a temporary
*   copy of each unit.
*/
inline unsigned Split::countOffBits() const
	{
    // This assumes that the unused bits in the last unit are all 0.
    unsigned num_bits_unset = 0;

    // Pretend that there are no excluded bits
    for (SplitTVect::const_iterator i = _unit.begin(); i != _unit.end(); ++i)
        {
        // Set comment in countOnBits about this Kernighan bit-counting algorithm
        split_t v = ~(*i);
        unsigned c = 0;
        for (; v; ++c)
            {
            v &= v - 1;
            }
        num_bits_unset += c;
        }

    if (!_excl_bits.empty())
        {
        // Correct count to account for excluded bits
        for (std::vector<unsigned>::const_iterator i = _excl_bits.begin(); i != _excl_bits.end(); ++i)
            {
            // Subtract one for every bit that is unset (remember that bits that are unset in this split
            // are set in its inverse)
            if (!isBitSet(*i))
                --num_bits_unset;
            }
        }

    // Account for irrelevant bits at the end
    unsigned num_irrelevant = (_nunits*_bits_per_unit - _split_ntax);
    num_bits_unset -= num_irrelevant;

    return num_bits_unset;
    }
	
/**
*   Returns the minimum of \c n (the number of taxa on one side of the split) and \c m (the number on the other side).
*   Trivial splits have \c m = 1 or \c n = 1, and thus have compexity 1, whereas the most complex split has complexity
*   #_split_ntax/2 (note that this maximum holds whether or not #_split_ntax is even or odd).
*/
inline unsigned Split::calcComplexity() const
	{
	const unsigned on_bits = countOnBits();
	const unsigned off_bits = countOffBits();
	return (on_bits < off_bits ? on_bits : off_bits);
	}
	
/**
*   Records the current values of the data members #_bits_per_unit, #_split_ntax and #_nunits, along with the size of a
*   split_t object on the current system, to a string object, which is returned.
*/
inline std::string Split::getDimensionInfo()
	{
    return boost::str(boost::format("Split information:\n  bits_per_unit = %d\n  split_ntax = %d\n  nunits = %d\n  sizeof(split_t) = %d")
        % _bits_per_unit % _split_ntax % _nunits % sizeof(split_t));
	}

/**
*   Returns true if and only if this object's #_unit is less than \c other's #_unit.
*/
inline bool Split::isLessThan(const Split & other) const
	{
    assert(_unit.size() == other._unit.size());
    return (_unit < other._unit);
	}

/**
*   Returns -1 if this object's #_unit is less than \c other's #_unit, 0 if this #_unit equals \c other's #_unit, and 1 if this #_unit
*   is greather than \c other's #_unit.
*/
inline int Split::cmp(const Split & other) const
	{
    assert(_unit.size() == other._unit.size());
    if (_unit < other._unit)
        return -1;
    else if (_unit == other._unit)
        return 0;
    else
        return 1;
	}

/**
*   Returns true if bit corresponding to taxon \c t is set, and returns false otherwise. 
*   The supplied value \c t should be 0 if it represents the first taxon.
*/
inline bool Split::isBitSet(unsigned t) const
	{
	assert(t < _split_ntax);
	
	unsigned i = t/_bits_per_unit;
	unsigned j = t % _bits_per_unit;
	split_t x = _unit[i] & (_split_unity << j);
	return (x > (split_t)0);
	}

/**
*   Deletes #_unit vector and creates a new #_unit data member with #_nunits elements all initialized to 0. Assumes
*   #_nunits is already set correctly, which will be the case if the member function calcNUnits() has just been
*   called.
*/
inline void Split::resize()
	{
    if (_unit.size() != _nunits)
        _unit.resize(_nunits, (split_t)0);
	}

/**
*   Sets the bit corresponding to taxon \c t, where \c t = 0 corresponds to the first taxon.
*/
inline void Split::setBit(unsigned t)
	{
	assert(t < _split_ntax);
	
	unsigned i = t/_bits_per_unit;
	unsigned j = t % _bits_per_unit;
	_unit[i] |= (_split_unity << j);
	}

/**
*   Sets all bits corresponding to values in the supplied vector \c bits_to_set. Values in \c bits_to_set are 0-based
*   indices (that is, use 0 to set the first bit).
*/
inline void Split::setBits(const std::vector<unsigned> bits_to_set)
	{
    for (std::vector<unsigned>::const_iterator it = bits_to_set.begin(); it != bits_to_set.end(); ++it)
        {
        unsigned value = *it;
	    assert(value < _split_ntax);
    	unsigned i = value/_bits_per_unit;
	    unsigned j = value % _bits_per_unit;
	    _unit[i] |= (_split_unity << j);
        }
	}

/**
*   Unsets the bit corresponding to taxon \c t, where \c t = 0 corresponds to the first taxon.
*/
inline void Split::unsetBit(unsigned t)
	{
	assert(t < _split_ntax);
	unsigned i = t/_bits_per_unit;
	unsigned j = t % _bits_per_unit;
	_unit[i] &= (~(_split_unity << j));
	}

/**
*   Unsets (that is, clears) all bits corresponding to values in the supplied vector \c bits_to_unset.
*   Values in \c bits_to_unset are 0-based indices (that is, use 0 to unset the first bit).
*/
inline void Split::unsetBits(const std::vector<unsigned> bits_to_unset)
	{
    for (std::vector<unsigned>::const_iterator it = bits_to_unset.begin(); it != bits_to_unset.end(); ++it)
        {
        unsigned value = *it;
	    assert(value < _split_ntax);
    	unsigned i = value/_bits_per_unit;
	    unsigned j = value % _bits_per_unit;
    	_unit[i] &= (~(_split_unity << j));
        }
	}

/**
*   If this split is subsumed in \c other, then a bitwise AND for any unit should equal the unit from this split. If this
*   tis not the case, it means there was a 0 in other for a bit that is set in this split. For example:
*
*	    **-*--**-*-  <-- this split
*	    *****-****-  <-- other split
*	    ----------------------------
*	    **-*--**-*-  <-- bitwise AND
*
*	In the above example, this split is subsumed in the other split because the bitwise AND equals this split. 
*	Note that every split is by definition subsumed in itself. The \c start_unit parameter is used by isCompatible().
*/
inline bool Split::subsumedIn(const Split & other, unsigned start_unit) const
	{
	for (unsigned i = start_unit; i < _nunits; ++i)
		{
		if ((_unit[i] & other._unit[i]) != _unit[i])
			return false;
		}
	return true;
	}

/**
*	Returns true if this split and \c other are compatible. The two splits a and b are compatible if a & b is nonzero
*   and also not equal to either a or b. For example, these two splits 
*
*       split a: -***---*--
*       split b: ----***--*
*         a & b: ----------
*
*   are compatible, because a & b = 0. The two splits below are also compatible because a & b == b:
*
*       split a: -****-*---
*       split b: --**--*--- <-
*         a & b: --**--*--- <-
*
*   These two splits, on the other hand, are not compatible because a & b != 0 and is not equal to either a or b:
*
*       split a: -***---*--
*       split b: ---***---*
*         a & b: ---*------
*/
inline bool Split::isCompatible(const Split & other) const
	{
	for (unsigned i = 0; i < _nunits; ++i)
		{
        split_t a       = _unit[i];
        split_t b       = other._unit[i];
		split_t a_and_b = (a & b);
        bool equals_a   = (a_and_b == a);
        bool equals_b   = (a_and_b == b);
        if (a_and_b && !(equals_a || equals_b))
            {
            // A failure of any unit to be compatible makes the entire split incompatible
            return false;
            }
		}

    // None of the units were incompatible, so that means the splits are compatible
	return true;
	}

/**
*	Converts each element of the _unit array to its bitwise complement, then clears the irrelevant bits in the final
*   element. The irrelevant bits arise because bits must be allocated in chunks of size #_bits_per_unit. For example, 
*   suppose #_split_ntax = 5, #_bits_per_unit = 4, and thus #_nunits = 2. The #_unit vector is shown below for a split
*   in which the bits for the last three taxa are "on":
*
*	    +---+---+---+---+  +---+---+---+---+
*	    | 0 | 0 | 0 | 1 |  | 1 | 1 | 0 | 0 |
*	    +---+---+---+---+  +---+---+---+---+
*	         _unit[1]          _unit[0]
*
*	The three most significant (that is, left-most) bits in _unit[1] in this case are irrelevant because there are 8 bits
*	total and only 5 of them are used to store information. Note that ordinarily #_bits_per_unit would be greater than 
*   4 (32 or 64 are common values for 32-bit or 64-bit systems, respectively). The #_unit vector shown above after 
*   invertSplit() is called would be:
*
*	    +---+---+---+---+  +---+---+---+---+
*	    | 0 | 0 | 0 | 0 |  | 0 | 0 | 1 | 1 |
*	    +---+---+---+---+  +---+---+---+---+
*	         _unit[1]          _unit[0]
*
*/
inline void Split::invertSplit()
	{
    for (SplitTVect::iterator i = _unit.begin(); i != _unit.end(); ++i)
		{
		split_t x = *i;
		*i = ~x;
		}

    // Unset the irrelevant bits at the end that do not correspond to any taxon
    split_t v = (split_t)(-1);              // v = 1111 (for example above); assumes split_t is an unsigned integer type
    v <<= (_split_ntax % _bits_per_unit);   // v = 1110 (introduce zeros for bits that are used)
    _unit[_nunits - 1] &= ~v;               // unit[1] = 1110, ~v = 0001, unit[1] & ~v = 0000

    // Must ensure that none of the bits now set corresponds to an excluded bit (cannot use the fast version of
    // countOnBits() unless all excluded bits have been cleared)
    if (!_excl_bits.empty())
        unsetBits(_excl_bits);
	}

/**
*	For each split in the \c ref_tree, checks to see if the split is in the \c test_tree SplitSet. If not, inserts the
*	split into \c missing. Assumes \c ref_tree and \c test_tree SplitSet objects are already filled.
*/
inline unsigned findSplitsAbsentInTestTree(const Split::SplitSet & ref_tree, const Split::SplitSet & test_tree, Split::SplitSet & missing)
	{
    missing.clear();
	Split::SplitSet::iterator insertLoc = missing.begin();
	unsigned num_absent_splits = 0;
	for (Split::SplitSet::const_iterator it = ref_tree.begin(); it != ref_tree.end(); ++it)
		{
		if (test_tree.find(*it) == test_tree.end())
			{
			insertLoc = missing.insert(insertLoc, *it);
			++num_absent_splits;
			}
		}
	return num_absent_splits;
	}

/**
*   Returns true if and only if this split is not equal to \c other.
*/
inline bool Split::operator!=(const Split & other) const
    {
    if (*this == other)
        return false;
    else
        return true;
    }

/**
*   Returns a split in which each element of the #_unit vector is its bitwise AND with the corresponding element from
*   the #_unit vector of \c other.
*   For example:
*
*       split a: -*-**-***
*       split b: -****----
*         a & b: -*-**----
*/
inline Split Split::operator&(const Split & other) const
	{
	Split tmp(*this);
	return tmp &= other;
	}
	
/**
*   Returns a reference to this Split object after replacing each element of its #_unit vector with the bitwise AND of
*   that element with the corresponding element of the unit vector of \c other. 
*   For example,
*
*       split a: -*-**-***
*       split b: -****----
*             a: -*-**---- after a &= b
*
*/
inline Split & Split::operator&=(const Split & other)
	{
	for (unsigned i = 0; i < _nunits; ++i)
		_unit[i] &= other._unit[i];
	return *this;
	}

/**
*   Returns a split in which each element of the #_unit vector is its bitwise XOR with the corresponding element from
*   the #_unit vector of \c other. 
*   For example:
*
*       split a: -*-**-***
*       split b: -****----
*         a ^ b: --*---***
*
*/
inline Split Split::operator^(const Split & other) const
    {
    Split tmp(*this);
    return tmp ^= other;
    }

/**
*   Returns a reference to this Split object after replacing each element of its #_unit vector with the bitwise XOR of
*   that element with the corresponding element of the unit vector of \c other.
*   For example:
*
*       split a: -*-**-***
*       split b: -****----
*             a: --*---*** after a ^= b
*
*/
inline Split & Split::operator^=(const Split & other)
    {
    for (unsigned i = 0; i < _nunits; ++i)
        _unit[i] ^= other._unit[i];
    return *this;
    }

/**
*   Returns a split in which each element of the #_unit vector is its bitwise OR with the corresponding element from
*   the #_unit vector of \c other.
*   For example:
*
*       split a: -*-**-***
*       split b: -****----
*         a | b: -****-***
*/
inline Split Split::operator|(const Split & other) const
	{
	Split tmp(*this);
	return tmp |= other;
	}

/**
*   Returns a reference to this Split object after replacing each element of its #_unit vector with the bitwise OR of
*   that element with the corresponding element of the unit vector of \c other.
*   For example:
*
*       split a: -*-**-***
*       split b: -****----
*             a: -****-*** after a |= b
*/
inline Split & Split::operator|=(const Split & other)
	{
	for (unsigned i = 0; i < _nunits; ++i)
		_unit[i] |= other._unit[i];
	return *this;
	}

/**
*	Calls isLessThan() for \c other and returns the value returned by that function.
*/
inline bool Split::operator<(const Split & other) const
	{
	return isLessThan(other);
	}

/**
*	Calls equals() on \c other and returns the value returned by that function.
*/
inline bool Split::operator==(const Split & other) const
	{
	return equals(other);
	}

/**
*	Returns true if all elements of the _unit vector are equal to all elements of the Split object \c other.
*/
inline bool Split::equals(const Split & other) const
	{
	return std::equal(_unit.begin(), _unit.end(), other._unit.begin());
	}

/**
*	Performs a bitwise OR of each element of _unit with the corresponding element in the _unit vector of \c other. The
*   result is a union of the sets defined by the two Split objects, and is useful in creating the Split for an interior
*   node, which is the union of the splits of its immediate descendants.
*/
inline void Split::combineWith(const Split & other)
	{
	*this |= other;
	}

/**
*	Converts this split into the intersection between this split and the other split. The other split remains unaffected
*	but this Split object is modified such that each element in the _unit array undergoes a bitwise AND operation with
*	the corresponding element from unit vector of \c other.
*/
inline void Split::intersectWith(const Split & other)
	{
	*this &= other;
	}

/**
*	Builds vector of unsigned integers each of which represents a bit that is "on" (that is, set to 1) in this split. The
*   value 0 refers to the first bit. This function removes values corresponding to excluded taxa before returning. This
*   function is utilized by both getOnList() and createNewickRepresentation(). Note that the supplied vector
*   \c v is cleared at the beginning of this function.
*/
inline void Split::getOnListImpl(std::vector<unsigned> & v) const
    {
    v.clear();
    unsigned k = 0;
    for (unsigned i = 0; i < _nunits; ++i)
        {
        for (unsigned j = 0; j < _bits_per_unit; ++j)
            {
            split_t bit = (_split_unity << j);
            bool is_on = ((_unit[i] & bit) > (split_t)0);
            if (is_on)
                v.push_back(k);
            if (++k == _split_ntax)
                break;
            }
        }

    // Eliminate values also in excl_bits vector
    if (!_excl_bits.empty())
        {
        std::vector<unsigned>::iterator last = std::set_difference(v.begin(), v.end(), _excl_bits.begin(), _excl_bits.end(), v.begin());
        v.erase(last, v.end());
        }
    }

/**
*	Returns vector of unsigned integers each of which represents a bit that is "on" (that is, set to 1) in this split. The
*   value 0 refers to the first bit. This function removes values corresponding to excluded taxa before returning the
*   vector.
*/
inline std::vector<unsigned> Split::getOnList() const
    {
    std::vector<unsigned> v;
    getOnListImpl(v);
    return v;
    }

/**
*	Builds vector of unsigned integers each of which represents a bit that is "off" (that is, set to 0) in this split. The
*   value 0 represents the first bit. This function removes values corresponding to excluded taxa before returning. Note
*   that the supplied vector \c v is cleared at the beginning of this function.
*/
inline void Split::getOffListImpl(std::vector<unsigned> & v) const
    {
    v.clear();
    unsigned k = 0;
    for (unsigned i = 0; i < _nunits; ++i)
        {
        for (unsigned j = 0; j < _bits_per_unit; ++j)
            {
            split_t bit = (_split_unity << j);
            bool is_on = ((_unit[i] & bit) > (split_t)0);
            if (!is_on)
                v.push_back(k);
            if (++k == _split_ntax)
                break;
            }
        }

    // Eliminate values also in excl_bits vector
    if (!_excl_bits.empty())
        {
        std::vector<unsigned>::iterator last = std::set_difference(v.begin(), v.end(), _excl_bits.begin(), _excl_bits.end(), v.begin());
        v.erase(last, v.end());
        }
    }

/**
*	Returns vector of unsigned integers each of which represents a bit that is "off" (that is, set to 0) in this split. The
*   value 0 represents the first bit. This function removes values corresponding to excluded taxa before returning the
*   list.
*/
inline std::vector<unsigned> Split::getOffList() const
    {
    std::vector<unsigned> v;
    getOffListImpl(v);
    return v;
    }

/**
*	Returns vector of unsigned integers each of which represents a bit that is currently excluded (first is 0). The bits
*   that are excluded are not considered either on or off, and are stored in the data member #_excl_bits.
*/
inline std::vector<unsigned> Split::getExcludedList() const
    {
    std::vector<unsigned> v;
    if (_excl_bits.empty())
        return v;
    else
        {
        v.resize(_excl_bits.size(), 0);
        std::copy(_excl_bits.begin(), _excl_bits.end(), v.begin());
        }
    return v;
    }

/**
*	Establishes a list of bits that are currently excluded. These bits will always be off, but they should not be 
*   returned in the list of values produced by getOffList(). Note that this function deletes any existing elements in the
*   #_excl_bits vector.
*/
inline void Split::setExcluded(const std::vector<unsigned> excl)
    {
    if (excl.empty())
        {
        _excl_bits.clear();
        }
    else
        {
        _excl_bits.resize(excl.size(), 0);
        std::copy(excl.begin(), excl.end(), _excl_bits.begin());

        // Must ensure _excl_bits is sorted because binary_search algorithm (used in createAndAppendPatternRepresentation())
        // and set_difference algorithm (used in getOnList() and getOffList()) require it
        std::sort(_excl_bits.begin(), _excl_bits.end());

        // Also ensure that none of the bits now set corresponds to an excluded bit (cannot use the fast version of
        // countOnBits() unless all excluded bits have been cleared)
        unsetBits(_excl_bits);
        }
    }

/**
*   Set the value of the data member #_on_symbol used to represent set bits by functions such as 
*   createPatternRepresentation().
*/
inline void Split::setOnSymbol(const char c)
    {
    _on_symbol = c;
    }

/**
*   Set the value of the data member #_off_symbol used to represent unset bits by functions such as 
*   createPatternRepresentation().
*/
inline void Split::setOffSymbol(const char c)
    {
    _off_symbol = c;
    }

/**
*   Set the value of the data member #_excl_symbol used to represent excluded bits by functions such as 
*   createPatternRepresentation().
*/
inline void Split::setExcludedSymbol(const char c)
    {
    _excl_symbol = c;
    }

/**
*   Returns current value of the data member #_on_symbol used to represent set bits by functions such as 
*   createPatternRepresentation().
*/
inline char Split::getOnSymbol() const
    {
    return _on_symbol;
    }

/**
*   Returns current value of the data member #_off_symbol used to represent unset bits by functions such as 
*   createPatternRepresentation().
*/
inline char Split::getOffSymbol() const
    {
    return _off_symbol;
    }

/**
*   Returns current value of the data member #_excl_symbol used to represent excluded bits by functions such as 
*   createPatternRepresentation().
*/
inline char Split::getExcludedSymbol() const
    {
    return _excl_symbol;
    }

/**
*   Returns current value of the data member #_split_ntax, which is the total number of taxa (counting both included and 
*   excluded bits) that can be represented by this Split object.
*/
inline unsigned Split::getNTaxa() const
    {
    return _split_ntax;
    }

/**
*   Sets the value of the data member #_split_ntax to the supplied value \c n, which results in the complete removal of
*   all information previously stored in the Split object (the member function clear() is called even if \c n is identical
*   to the current value of #_split_ntax).
*/
inline void Split::setNTaxa(unsigned n)
    {
    clear();
    calcNUnits(n);
    resize();
    }

}   //namespace strom

#endif //defined(GALAX_SPLIT_HPP)

