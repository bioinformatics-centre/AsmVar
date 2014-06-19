map<size_t,size_t> AxtVar::FindSeed ( vector<MapReg> & mapreg ) {

    map<int, map<size_t,size_t> > mark;
    map<size_t,size_t> tmpmap;
    if ( mapreg.size() > 0 ) {

        tmpmap[0]=0;
        unsigned int prepos( mapreg[0].target.start );
        for ( size_t i(1); i < mapreg.size(); ++i ) {

            if ( mapreg[i].target.start < prepos ) {
                mark.insert( make_pair(tmpmap.size(), tmpmap) ); tmpmap.clear(); // For 'mark' if key conflict take the last one
            }
            tmpmap[i] = i;
            prepos    = mapreg[i].target.start;
        }
        mark.insert( make_pair(tmpmap.size(), tmpmap) ); tmpmap.clear();

        size_t j = mapreg.size() - 1; // The last element
        prepos   = mapreg[j].target.start;
        tmpmap[j]= j;
        for ( --j ; j >= 0; --j ) {
            if ( mapreg[j].target.start < prepos ) {
                mark.insert( make_pair(tmpmap.size(), tmpmap) ); tmpmap.clear(); // For 'mark' if key conflict take the last one
            }
            tmpmap[j] = j;
            prepos    = tmpmap[j].target.start;
        }
        mark.insert( make_pair(tmpmap.size(), tmpmap) ); tmpmap.clear();
        tmpmap = mark.rbegin()->second();
    }

    return tmpmap;
}
