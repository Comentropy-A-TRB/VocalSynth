#include<bits/stdc++.h>
using namespace std;

const int tabDelta[12]={0,2,-9,-7,-5,-4,-2,0};
const int A1=33;    
const double unit=327.868;   // MIDI 文件中的 tempo 指的是每拍微妙 (microseconds) 数。这里乘一个 1/1000
int mov=1;
int keyToNote(const string &s){    // A1 -> 55Hz
    // cerr<<s<<" "<<s.size()<<endl;
    int up=0;
    if(s[0]=='#')   up=1;
    assert((int)s.size()==2+up);
    assert((int)s.size()==2+up && isalpha(s[up]) && isdigit(s[up+1]));
    return A1+tabDelta[s[up]-'A']+up+12*(s[up+1]-'1')+12*mov;
}

string s[]={"C","#C","D","#D","E","F","#F","G","#G","A","#A","B"};
string t[]={"1","2","3","4","5","6","7","8"};

char lyr[2000];
map<int,string> tab;

void print(const string &hit,const string &str,double beat){
    printf("%s %.6lf %s\n0\n0\n0\n",hit.c_str(),beat,str.c_str());
}
bool is_vow(char c){
    return c=='a' || c=='i' || c=='u' || c=='e' || c=='o' || c=='n';
}
void split_lyric(string &hit,double beat){   
    string vow="";
    int i=strlen(lyr)-1;
    bool flg=0;
    for(;i>=0 && is_vow(lyr[i]);i--){
        if(lyr[i]=='n' && flg)  break ;
        vow+=lyr[i],flg|=(lyr[i]!='n');
    }
    reverse(vow.begin(),vow.end());
    if(~i){
        string con="";
        double len=min(beat*0.5,0.1);
        for(;i>=0;i--)
            con+=lyr[i];
        reverse(con.begin(),con.end());
        print(hit,con,len);
        beat-=len;
    }
    print(hit,vow,beat);
}
int main(){
    // ofstream debug("debug.txt");
    freopen("res.txt","r",stdin);
    freopen("reyichang.txt","w",stdout);

    for(int i=0;i<12;i++)
        for(int j=0;j<8;j++){
            string now=s[i]+t[j];
            // cerr<<i<<" "<<j<<" "<<now<<endl;
            int res=keyToNote(now);
            if(0<=res && res<128)
                tab[res]=now;
        }
    string tmp;
    int nam=0,tim=0;
    string hit;
    auto readMeta=[&](){
        int pos=tmp.find("text=")+6,idx=0;
        while(tmp[pos]!='\'')   lyr[idx++]=tmp[pos++];
        lyr[idx]=0;
        tim=0;
        for(int i=pos;i<tmp.size();i++)
            if(isdigit(tmp[i])) tim=tim*10+tmp[i]-'0';
        return ;
    };
    while(getline(cin,tmp)){
        if(tmp.find("lyrics")!=string::npos){
            // 这么写……只是读到空白符再终止的……我要无语了。
            // sscanf(tmp.c_str(),"MetaMessage(\'lyrics\', text=\'%s\', time=%d)",lyr,&tim);
            // lyr[strlen(lyr)-2]=0;
            readMeta();
            cerr<<lyr<<" "<<tim<<"\n";
            if(tim){
                print("A3","SP",tim/480.0); // ticks 数。
            }
        }else if(tmp.find("note_on")!=string::npos){
            sscanf(tmp.c_str(),"note_on channel=%*d note=%d velocity=%*d time=%*d",&nam);
            hit=tab[nam];
        }else if(tmp.find("note_off")!=string::npos){
            sscanf(tmp.c_str(),"note_off channel=%*d note=%*d velocity=%*d time=%d",&tim);
            split_lyric(hit,tim/480.0);
        }
    }
    return 0;
}