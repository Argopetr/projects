#ifndef CLISTL2_H
#define CLISTL2_H

#include <iostream>
#include "CListBuffer.h"

#define START_LEN 2
#define COMBO 2

#define LIST_RANGE_ERROR -1

template <typename CItem>
struct CListL2 {
private:
    int *OrderArr, Amount, Size;
    struct CListCluster {
    private:
        int Size, Amount;
        CItem *Space;

    public:
        CListCluster *Back, *Next;

        CListCluster();
        CListCluster(int size);
        CListCluster(const CListCluster& cls);
        ~CListCluster() {delete[] Space;}

        int del(int ind);
        CItem add(const CItem& val);
        int search(const CItem& val) const;

        int size() const {return Size;}
        int amount() const {return Amount;}
        int isEmpty(void) const {return (Amount)? 0 : 1;}

        CItem getVal(int ind) const;
    } *Source, *Head;

public:
    CListL2();
    CListL2(const CListL2& lst);
    ~CListL2() {clean();}

    CListL2& add(const CItem& val);
    CListL2& add(int ind, const CItem& val);
    CListL2& del(int ind);

    int isEmpty(void) const {return (Amount)? 0: 1;}
    int amount(void) const {return Amount;}

    CItem getVal(int ind) const;
    CListL2& setVal(int ind, const CItem& val);

    CListL2& join(const CListL2& lst);
    CListL2& copy(const CListL2& lst);
    CListL2& clean(void);

    int search(const CItem& val) const;
    CListL2& sortInc(void);
    CListL2& sortDec(void);

    CListBuffer<CListL2, CItem> operator[](int ind);
    CListL2& operator=(const CListL2& lst) {return copy(lst);}
    CListL2& operator+=(const CListL2& lst) {return join(lst);}
    CListL2& operator+=(const CItem& val) {return add(val);}
    CListL2 operator+(const CListL2& lst) const;
    CListL2 operator+(const CItem& val) const;
};

// CLUSTER OF LIST

template <typename CItem>
CListL2<CItem>::CListCluster::CListCluster() {
    Size = 0;
    Amount = 0;

    Space = NULL;

    Back = NULL;
    Next = NULL;
}

template <typename CItem>
CListL2<CItem>::CListCluster::CListCluster(int size) {
    if ( size > 0 ) {
        Size = size;
        Space = new CItem[size];
    } else {
        Size = 0;
        Space = NULL;
    }

    Amount = 0;

    Back = NULL;
    Next = NULL;
}

template <typename CItem>
CListL2<CItem>::CListCluster::CListCluster(const CListCluster &cls) {
    Size = cls.Size;
    Amount = cls.Amount;

    Space = new CItem[Size];

    for(int i = 0; i < Amount; i++)
        Space[i] = cls.Space[i];

    Back = NULL;
    Next = NULL;
}

template <typename CItem>
int CListL2<CItem>::CListCluster::del(int ind) {
    int err = 0;

    if ( ind >= 0 && ind < Amount ) {
        Amount--;

        for(int i = ind; i < Amount; i++)
            Space[i] = Space[i + 1];
    } else
        err = -1;

    return err;
}

template <typename CItem>
CItem CListL2<CItem>::CListCluster::add(const CItem &val) {
    CItem res = val;

    if ( Size ) {
        if ( Amount ) {
            if ( Amount < Size ) {
                if ( val >= Space[Amount - 1] )
                    Space[Amount] = val;
                else {
                    int ind = 0, delta = 1;
                    CItem tmp;

                    for(; delta < Amount - 1; delta <<= 1);

                    for(; delta > 0; delta >>= 1)
                        if ( ind + delta < Amount && val >= Space[ind + delta] )
                            ind += delta;

                    if ( val >= Space[ind] )
                        ind++;

                    for(int i = ind; i <= Amount; i++) {
                        tmp = Space[i];
                        Space[i] = res;
                        res = tmp;
                    }
                }

                Amount++;
            } else if ( val < Space[Amount - 1] ) {
                int ind = 0, delta = 1;
                CItem tmp;

                for(; delta < Amount - 1; delta <<= 1);

                for(; delta > 0; delta >>= 1)
                    if ( ind + delta < Amount && val >= Space[ind + delta] )
                        ind += delta;

                if ( val >= Space[ind] )
                    ind++;

                for(int i = ind; i < Amount; i++) {
                    tmp = Space[i];
                    Space[i] = res;
                    res = tmp;
                }
            }
        } else {
            Space[0] = val;
            Amount = 1;
        }
    }

    return res;
}

template <typename CItem>
int CListL2<CItem>::CListCluster::search(const CItem &val) const {
    int res = -1;

    if ( Amount && val >= Space[0] && val <= Space[Amount - 1] ) {
        int ind = 0, delta = 1;

        for(; delta < Amount - 1; delta <<= 1);

        for(; delta > 0; delta >>= 1)
            if ( ind + delta < Amount && val >= Space[ind + delta] )
                ind += delta;

        if ( val == Space[ind] )
            res = ind;
    }

    return res;
}

template <typename CItem>
CItem CListL2<CItem>::CListCluster::getVal(int ind) const {
    CItem res;

    if ( ind >= 0 && ind < Amount )
        res = Space[ind];
    else {
        std::cerr << "Out of list memory" << std::endl;
        exit(LIST_RANGE_ERROR);
    }

    return res;
}

// LIST

template <typename CItem>
CListL2<CItem>::CListL2() {
    Size = 0;
    Amount = 0;

    OrderArr = NULL;
    Source = NULL;
    Head = NULL;
}

template <typename CItem>
CListL2<CItem>::CListL2(const CListL2& lst) {
    Amount = lst.Amount;
    Size = lst.Size;
    OrderArr = new int[Size];

    for(int i = 0; i < Amount; i++)
        OrderArr[i] = lst.OrderArr[i];

    if ( Amount ) {
        CListCluster *tmp = lst.Source, *bff;

        Source = new CListCluster(*lst.Source);
        Head = Source;

        for(; tmp -> Next ;) {
            bff = new CListCluster(*tmp -> Next);

            Head -> Next = bff;
            Head = bff;

            tmp = tmp -> Next;
        }
    } else {
        Source = NULL;
        Head = NULL;
    }
}

template <typename CItem>
CListL2<CItem>& CListL2<CItem>::add(const CItem &val) {
    if ( Amount ) {
        CListCluster *tmp = Source;
        CItem tmpVal = val;
        int num, startInd = 0, flag;

        if ( Amount == Size ) {
            CListCluster *newHead = new CListCluster(Size * (COMBO - 1) + START_LEN);
            int *newOrderArr = new int[Size * COMBO + START_LEN];


            Head -> Next = newHead;
            newHead -> Back = Head;
            Head = newHead;

            for(int i = 0; i < Size; i++)
                newOrderArr[i] = OrderArr[i];

            delete[] OrderArr;
            OrderArr = newOrderArr;

            Size = Size * COMBO + START_LEN;
        }

        for(; tmp ;) {
            tmpVal = tmp -> add(tmpVal);
            flag = tmp -> search(val);

            if ( flag != -1 )
                num = startInd + flag;

            startInd = startInd * COMBO + START_LEN;
            tmp = tmp -> Next;
        }

        for(int i = 0; i < Amount; i++)
            if ( OrderArr[i] >= num )
                OrderArr[i]++;

        OrderArr[Amount] = num;
        Amount++;
    } else {
        Size = START_LEN;
        Amount++;

        Source = new CListCluster(START_LEN);
        Head = Source;

        OrderArr = new int[START_LEN];

        Source -> add(val);
        OrderArr[0] = 0;
    }

    return *this;
}

template <typename CItem>
CListL2<CItem>& CListL2<CItem>::add(int ind, const CItem &val) {
    if ( Amount && ind >= 0 && ind <= Amount ) {
        CListCluster *tmp = Source;
        CItem tmpVal = val;
        int num, startInd = 0, flag;

        if ( Amount == Size ) {
            CListCluster *newHead = new CListCluster(Size * (COMBO - 1) + START_LEN);
            int *newOrderArr = new int[Size * COMBO + START_LEN];

            Head -> Next = newHead;
            newHead -> Back = Head;
            Head = newHead;

            for(int i = 0; i < Size; i++)
                newOrderArr[i] = OrderArr[i];

            delete[] OrderArr;
            OrderArr = newOrderArr;

            Size = Size * COMBO + START_LEN;
        }

        for(; tmp ;) {
            tmpVal = tmp -> add(tmpVal);
            flag = tmp -> search(val);

            if ( flag != -1 )
                num = startInd + flag;

            startInd = startInd * COMBO + START_LEN;
            tmp = tmp -> Next;
        }

        for(int i = Amount; i >= 0; i--) {
            if ( i > ind )
                OrderArr[i] = OrderArr[i - 1];

            if ( OrderArr[i] >= num )
                OrderArr[i]++;
        }

        OrderArr[ind] = num;
        Amount++;
    } else if ( !ind ) {
        Size = START_LEN;
        Amount = 1;

        Source = new CListCluster(START_LEN);
        Head = Source;

        OrderArr = new int[START_LEN];

        Source[0] = val;
        OrderArr[0] = 0;
    }

    return *this;
}

template <typename CItem>
CListL2<CItem>& CListL2<CItem>::del(int ind) {
    if ( ind >= 0 && ind < Amount ) {
        CListCluster *block = Source;
        int blockLen = START_LEN, num = OrderArr[ind];
        
        Amount--;

        for(int i = 0; i < Amount; i++) {
            if ( i >= ind )
                OrderArr[i] = OrderArr[i + 1];

            if ( OrderArr[i] > num )
                OrderArr[i]--;
        }

        for(; num >= blockLen ;) {
            num -= blockLen;
            blockLen = blockLen * COMBO;
            block = block -> Next;
        }

        for(; block ;) {
            block -> del(num);

            if ( block -> Next )
                block -> add(block -> Next -> getVal(0));

            block = block -> Next;
            num = 0;
        }
    }

    return *this;
}

template <typename CItem>
CItem CListL2<CItem>::getVal(int ind) const {
    CItem res;

    if ( ind >= 0 && ind < Amount ) {
        CListCluster *block = Source;
        int blockLen = START_LEN;
        int num = OrderArr[ind];

        for(; num >= blockLen ;) {
            num -= blockLen;
            blockLen = blockLen * COMBO;
            block = block -> Next;
        }

        res = block -> getVal(num);
    } else {
        std::cerr << "Out of list memory" << std::endl;
        exit(LIST_RANGE_ERROR);
    }

    return res;
}

template <typename CItem>
CListL2<CItem>& CListL2<CItem>::setVal(int ind, const CItem& val) {
    if ( ind >= 0 && ind < Amount ) {
        del(ind);
        add(ind, val);
    } else {
        std::cerr << "Out of list memory" << std::endl;
        exit(LIST_RANGE_ERROR);
    }

    return *this;
}

template <typename CItem>
CListL2<CItem>& CListL2<CItem>::clean(void) {
    CListCluster *tmp;

    for(; Source ;) {
        tmp = Source -> Next;

        delete Source;

        Source = tmp;
    }

    delete[] OrderArr;

    Size = 0;
    Amount = 0;

    Source = NULL;
    Head = NULL;
    OrderArr = NULL;

    return *this;
}

template <typename CItem>
CListL2<CItem>& CListL2<CItem>::copy(const CListL2& lst) {
    if ( &lst != this ) {
        clean();

        Amount = lst.Amount;
        Size = lst.Size;
        OrderArr = new int[Size];

        for(int i = 0; i < Amount; i++)
            OrderArr[i] =lst.OrderArr[i];

        if ( Amount ) {
            CListCluster *tmp = lst.Source, *bff;

            Source = new CListCluster(*lst.Source);
            Head = Source;

            for(; tmp -> Next ;) {
                bff = new CListCluster(*tmp -> Next);

                Head -> Next = bff;
                Head = bff;

                tmp = tmp -> Next;
            }
        } else {
            Source = NULL;
            Head = NULL;
        }
    }

    return *this;
}

template <typename CItem>
CListL2<CItem>& CListL2<CItem>::join(const CListL2& lst) {
    for(int i = 0; i < lst.Amount; i++)
        add(lst.getVal(i));


    return *this;
}

template <typename CItem>
int CListL2<CItem>::search(const CItem& val) const {
    CListCluster *block = Source;
    int res = -1, startBlockInd = 0;

    for(; block ;) {
        res = block -> search(val);

        if ( res != -1 ) {
            res += startBlockInd;
            block = NULL;
        } else {
            startBlockInd = startBlockInd * COMBO + START_LEN;
            block = block -> Next;
        }
    }

    if ( res != -1 )
        for(int i = 0; i < Amount; i++)
            if ( OrderArr[i] == res ) {
                res = i;
                i = Amount;
            }

    return res;
}

template <typename CItem>
CListL2<CItem>& CListL2<CItem>::sortInc(void) {
    for(int i = 0; i < Amount; i++)
        OrderArr[i] = i;

    return *this;
}

template <typename CItem>
CListL2<CItem>& CListL2<CItem>::sortDec(void) {
    for(int i = 0; i < Amount; i++)
        OrderArr[i] = Amount - i - 1;

    return *this;
}

template <typename CItem>
CListBuffer<CListL2<CItem>, CItem> CListL2<CItem>::operator[](int ind) {
    CListBuffer<CListL2<CItem>, CItem> res(this, ind);

    return res;
}

template <typename CItem>
std::ostream& operator<<(std::ostream& output, const CListL2<CItem>& lst) {
    if ( lst.amount() ) {
        output <<'[';

        for(int i = 0; i < lst.amount(); i++) {
            if ( i )
                output << ", ";

            output << lst.getVal(i);
        }

        output <<']';
    } else
        output << "Empty list";

    return output;
}

template <typename CItem>
std::istream& operator>>(std::istream& input, CListL2<CItem>& lst) {
    int num;
    CItem tmp;

    input >> num;

    for(; num > 0; num--) {
        input >> tmp;
        lst.add(tmp);
    }

    return input;
}

template <typename CItem>
CListL2<CItem> CListL2<CItem>::operator+(const CListL2<CItem>& lst) const {
    CListL2 res(*this);

    return res.join(lst);
}

template <typename CItem>
CListL2<CItem> CListL2<CItem>::operator+(const CItem& val) const {
    CListL2 res(*this);

    return res.add(val);
}

#undef START_LEN
#undef COMBO

#endif
