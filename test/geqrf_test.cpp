#include <vector>
#include <algorithm>

#include "ampclapack.h"

int main()
{
    int info;
    int nb = 128;
    int n = nb*10;

    std::vector<float> a(n*n, 1);

    // scale diagonal
    // for (int i = 0; i < (n*n); i += (n+1))
    //     a[i] *= n;

    amplapack_spotrf('U', n, &a.at(0), n, &info);

    printf("last = %f\n", a.at((n*n)-1));
}
