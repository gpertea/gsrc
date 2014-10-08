//patternInit:

_myersAdjustBitmask(PatternState_<xxx> &state, TValue const value, TShift const shift, False) {
    register unsigned ord = ordValue(value); //ascii code of character, for regular text
	register unsigned short x = shift - state.shift[ord];
	if (x < BitsPerValue<TWord>::VALUE)
		state.bitMasks[ord] = (state.bitMasks[ord] >> x) | ((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
	else
		state.bitMasks[ord] = (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
	state.shift[ord] = shift;
}

_myersGetBitmask(PatternState_<xxx> &state, TValue const value, TShift const shift, False) 
{
	typedef typename PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord TWord;

    if (IsSameType<TFinderCSP, NMatchesNone_>::VALUE && value == unknownValue<TValue>())
        return 0;

    if (IsSameType<TFinderCSP, NMatchesAll_>::VALUE && value == unknownValue<TValue>())
        return (shift < BitsPerValue<TWord>::VALUE)? -1 << shift: -1;
    
    register unsigned ord = ordValue(value);
    register TWord res;
    register TShift x = shift - state.shift[ord];
	if (x < BitsPerValue<TWord>::VALUE) 
		res = state.bitMasks[ord] >> x;
	else
		res = 0;

    if (IsSameType<TPatternCSP, NMatchesAll_>::VALUE)
    {
        ord = ordValue(unknownValue<TValue>());
        x = shift - state.shift[ord];
        if (x < BitsPerValue<TWord>::VALUE) 
            res |= state.bitMasks[ord] >> x;
    }
    return res;
}


//_patternInitSmallStateBanded
template <typename TFinder, typename TNeedle, typename TNeedle2, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
inline bool 
 _patternInitSmallStateBanded(
    TFinder &finder,
	TNeedle2 const & needle, 
    PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
	typedef Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TPattern;
	typedef typename TPattern::TWord TWord;
	typedef typename Iterator<TNeedle2 const, Standard>::Type TIter;
	typedef typename Value<TNeedle>::Type TValue;

	//_myersPreInit(state, typename MyersSmallAlphabet_<TValue>::Type());
    {
    memset(state.bitMasks, 0, (ValueSize<TValue>::VALUE + 1) * sizeof(state.bitMasks[0]));
    memset(state.shift, 0, ValueSize<TValue>::VALUE * sizeof(state.shift[0]));
    
    }

	typename Size<TNeedle>::Type const ndlLength = length(needle);
	
	// Initialize row 0 either with zeros or increasing numbers
	// This can be realized using the following DP pattern and
	// assuming character mismatches at rows -1, -2,...
	// Thus we initialize the bitmasks and VN with 0.
	// VP depends on global/local alignment
	//
	//  0  1  2  3  4   -2 -2 -2 -2 -2   (row -2)
	//  0  1  2  3  4   -1 -1 -1 -1 -1   (row -1)
	//  0  1  2  3  4    0  0  0  0  0   (row  0)
	//  1                1
	//	    global           local
	//
	//  VP = 100...      VP = 111...
	//

    //                     HP[0]==1 <-> global, HP[0]==0 <-> local
    register TWord VP = (MyersUkkonenHP0_<TSpec>::VALUE == 1)? (TWord)1 << ((int)BitsPerValue<TWord>::VALUE-1): maxValue<TWord>(); 
    register TWord VN = 0;
	
	// Errors are counted along the lowest diagonal and the
	// lowest row of the band.
	// The errors in the top-left corner are 0.
	//
	// 0 * * * *
	//   x * * * *
	//     x * * * *
	//       x x x x x
	//
	//       |-------|
	//     diagWidth + 1 = 5
	//
	// diagWidth = length(container(finder)) + state.leftClip + state.rightClip - length(needle)
	
    register unsigned errors = 0;
	TIter ndlIter = begin(needle, Standard()); //pointer at the beginning of needle
	TIter ndlEnd;
	
	// The errors along the diagonal can only increase or stay the same.
	// There is only the last row of length diagWidth where errors can decrease.
	// If errors exceeds cutOff it cannot reach maxErrors again.
	
	
	typename Size<TFinder>::Type const columns = length(container(hostIterator(finder))) + state.leftClip;
    //columns = strlen(T)
    register unsigned cutOff = state.maxErrors; // = maximum edit distance
	if (columns > ndlLength)
	{
		cutOff += columns - ndlLength;		// clipping case *0
		ndlEnd = end(needle, Standard()); //ndlEnd = pointer to the '\0' at the end of needle
	} else { //pattern shorter than text
        errors += ndlLength - columns;
		ndlEnd = ndlIter + columns;			// clipping case *1
	}

	unsigned short shift = 0;
	
	if (state.leftClip != 0)
	{
		//////////////////////////////////////////////////////////////////
		// PART 0: go down the parallelogram in a empty (clipped) area
		//////////////////////////////////////////////////////////////////

		errors += state.leftClip;
		if (errors > ndlLength) errors = ndlLength;
		if (errors > cutOff) return false;

	// leftClip = 2
	//   |-|
	//
	//   . . * * *
	//     . * * * *
	//       * * * * *
	//         * * * * .
	//           * * * . .
	//             * * . . .
	//
	//                 |---|
	//               rightClip = 3
    //
	// We divide the parallelogam into 3 sections:
	//
	//   A A A A
	//     A A A B
	//       A A B B
	//         A B B C
	//           B B C C
	//             B C C C
	//               C C C C
	//
	// Depending on where the clipping ends we identify 4 different clipping cases:
	//
	//	 case 00            case 10            case 01            case 11
	//   . . * *            . . . .            . . * *            . . . .
	//     . * * *            . . . *            . * * *            . . . *
	//       * * * *            . . * *            * * * .            . . * .
	//         * * * *            . * * *            * * . .            . * . .
	//           * * * .            * * * .            * . . .            * . . .
	//             * * . .            * * . .            . . . .            . . . .
	//

		// adjust bitmasks (errors = number of needle chars to preprocess)
		for (; shift < errors; ++ndlIter, ++shift)
			_myersAdjustBitmask(state, getValue(ndlIter), shift, typename MyersSmallAlphabet_<TValue>::Type());

		// initialise left column with
		//
		//  0  1  2  3  4   -2 -2 -2 -2 -2
		//  0  1  2  3  4   -1 -1 -1 -1 -1
		//  0  1  2  3  4    0  0  0  0  0
		//  1                1
		//  2   global       2   local
		//  3                3
		//  4                4
		//
		//  VP = 111100...   VP = 111111...
		if (errors < (unsigned)BitsPerValue<TWord>::VALUE-1)
			VP |= ((TWord) -1) << ((unsigned)BitsPerValue<TWord>::VALUE-1 - errors);
		else
			VP = (TWord)-1;
	}

	for (; ndlIter != ndlEnd; ++ndlIter, goNext(finder), ++shift)
	{
		//////////////////////////////////////////////////////////////////
		// PART 1: go down the parallelogram
		//////////////////////////////////////////////////////////////////
	
	// adjust bitmask
	_myersAdjustBitmask(state, getValue(ndlIter), shift, typename MyersSmallAlphabet_<TValue>::Type());
	
		/////////////////////////
	// DIAGONAL MYERS CORE
		
		// VP/VN --> D0  (original Myers)
	    register TWord X = _myersGetBitmask(state, ordValue(*finder), shift, typename MyersSmallAlphabet_<TValue>::Type()) | VN;
	    register TWord D0 = ((VP + (X & VP)) ^ VP) | X;
		
		// adjust errors corresponding to rightmost bit of D0
        errors += (~D0 >> (BitsPerValue<TWord>::VALUE - 1)) & 1;
        if (errors > cutOff) return false;

		// D0 --> HP/HN  (original Myers)
        register TWord HN = VP & D0;
        register TWord HP = VN | ~(VP | D0);
    //    const int PADDING = sizeof(TWord)*2 + 1;
    //    std::cerr << std::hex;
    //    std::cerr << "\tD0"<<std::setw(PADDING)<<(__uint64)D0<<"\tHN"<<std::setw(PADDING)<<(__uint64)HN<<"\tHP"<<std::setw(PADDING)<<(__uint64)HP << std::endl;

		// moving register down corresponds to shifting HP/HN up (right shift)
		// HP/HN --> shift --> VP/VN (modified Myers)
        X = D0 >> 1;
        VN = X & HP;
        VP = HN | ~(X | HP);
    //    std::cerr << "\t  "<<std::setw(PADDING)<<' '<<"\tVN"<<std::setw(PADDING)<<(__uint64)VN<<"\tVP"<<std::setw(PADDING)<<(__uint64)VP << std::endl;
    //    std::cerr << std::dec;
    
    } //for
    state.VP0 = VP;
    state.VN0 = VN;
    state.errors = errors;
    _myersPostInit(state, typename MyersSmallAlphabet_<TValue>::Type());
    return true;
}



//inline bool _findMyersSmallPatternsBanded(
	TFinder & finder, 
	TNeedle const & needle,
    PatternState_<TNeedle2, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
	typedef PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TState;
	typedef typename TState::TWord TWord;
	typedef typename Iterator<TNeedle const, Standard>::Type TIter;
	typedef typename Value<TNeedle>::Type TValue;

    register TWord VP = state.VP0;
    register TWord VN = state.VN0;
    register TWord errors = state.errors;
    register TWord const maxErrors = state.maxErrors;
	register unsigned short const shift = length(needle);

	for (; !atEnd(finder); goNext(finder))
    {
		// PART 2: go right
		// normal Myers
		register TWord X = _myersGetBitmask(state, ordValue(*finder), shift, typename MyersSmallAlphabet_<TValue>::Type()) | VN;
		register TWord D0 = ((VP + (X & VP)) ^ VP) | X;
		register TWord HN = VP & D0;
		register TWord HP = VN | ~(VP | D0);
		X = (HP << 1) | 1;
		VN = X & D0;
		VP = (HN << 1) | ~(X | D0);
        errors += (HP >> (BitsPerValue<TWord>::VALUE - 2)) & 1;
        errors -= (HN >> (BitsPerValue<TWord>::VALUE - 2)) & 1;
                if (errors <= maxErrors)
        {
            state.VP0 = VP;
            state.VN0 = VN;
            state.errors = errors;
			_setFinderEnd(finder);
			if (IsSameType<TSpec, FindPrefix>::VALUE)
			{
				_setFinderLength(finder, endPosition(finder));
			}
            return true;
        }
    }
