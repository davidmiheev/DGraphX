#include "Interp2D.h"
#include "DGraphX.h"

void addpoint(char s[], int p);
void reverse(char s[]);

void addpoint(char s[], int p) {
    //char dest[1]; dest[0] = '0';
    int j = strlen(s);
    for(int i = strlen(s); i > j - p - 1; i--) s[i + 1] = s[i];
    s[j - p] = '.'; j = strlen(s);
    if(s[0] == '.') { for(int i = strlen(s); i > j - p - 2; i--) s[i + 1] = s[i];
        s[0] = '0';
    }
    if(s[1] == '.' && s[0] == '-') { for(int i = strlen(s); i > j - p - 2; i--) s[i + 1] = s[i];
        s[1] = '0';
    }
}
void reverse(char s[])
{
    int i, j; char c;
    
    for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}

void dtoa(double x, char s[]) {
    int i, sign, n = (int) 100*x;
    if ((sign = n) < 0)
        n = -n; i = 0;
    do {
        s[i++] = n % 10 + '0';
    } while ((n /= 10) > 0);
    if (sign < 0)
        s[i++] = '-';
    s[i] = '\0';
    reverse(s);
    addpoint(s,2);
}

void itoa(int n, char s[]) {
    int i, sign;
    
    if ((sign = n) < 0)
        n = -n;
    i = 0;
    do {
        s[i++] = n % 10 + '0';   
    } while ((n /= 10) > 0);
    if (sign < 0)
        s[i++] = '-';
    s[i] = '\0';
    reverse(s);
}

void scat(char* s, char* str, unsigned long len, char* msg) {

    for (unsigned long i = 0; i < len + strlen(str); i++) {
        if(i < len) msg[i] = s[i];
        else msg[i] = str[i - len];
    } msg[len + strlen(str)] = '\0';
}
/*int main() {
    char s[128];
    dtoa(25.08457,s);
    printf("%s\n",s);
    return 0;
}*/
