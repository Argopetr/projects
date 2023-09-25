#ifndef CLISTBUFFER_H
#define CLISTBUFFER_H

template <typename CList, typename CItem>
struct CListBuffer {
    int Ind;
    CList *List;

    CListBuffer(CList *lst, int ind);

    operator CItem() {return List -> getVal(Ind);}
    CItem operator=(const CItem& val);
};

template <typename CList, typename CItem>
CListBuffer<CList, CItem>::CListBuffer(CList *lst, int ind) {
    List = lst;
    Ind = ind;
}

template <typename CList, typename CItem>
CItem CListBuffer<CList, CItem>::operator=(const CItem& val) {
    List -> setVal(Ind, val);

    return val;
}

#endif
