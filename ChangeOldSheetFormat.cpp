#include<bits/stdc++.h>
using namespace std;

int main(){
    ifstream she("daisy.txt");
    ofstream now("new.txt");
    string s;
    while(getline(she,s)){
        now<<s<<endl;
        
        for(int i=0;i<3;i++)
            now<<"0"<<endl;
        now<<endl;
    }
    return 0;
}