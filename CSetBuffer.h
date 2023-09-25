#ifndef CSETBUFFER_H
#define CSETBUFFER_H

template <typename CSet, typename CItem>
struct CSetBuffer {
private:
	CSet *Set;
    const CItem *Item;
    int Flip;

public:
    CSetBuffer(CSet* set, const CItem& item);

	operator int() {return Flip;}
    int operator=(int val);
};

template <typename CSet, typename CItem>
CSetBuffer<CSet, CItem>::CSetBuffer(CSet* set, const CItem& item) {
    Set = set;
    Item = &item;
    Flip = set -> contains(item);
}

template <typename CSet, typename CItem>
int CSetBuffer<CSet, CItem>::operator=(int val) {
    if (val)
        Set -> add(*Item);
    else
        Set -> del(*Item);

    return Flip = val;
}

#endif
