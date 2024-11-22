# VocalSynth

This is an attempt to synthesize vocal without sampling.

Welcome everybody to improve this algorithm.

This reponsitrory is based on ```FFTW``` and [AudioFile](https://github.com/adamstark/AudioFile).

## Progress

Now, the author has finished the function of synthesizing vowels and some consonants.

## Wish

Want to synth consonants like 'l','r','n' better.

## Basic theory

### Vowels

懒得写英文了，简单来说就是，一个音高上的音，由基频和泛音列组成，基频（f0）表征音高，而泛音列表征音色，泛音的频率会是基频的正整数倍，这可以参照高中物理。元音的合成基于如下事实：由于某些特定发声部位的共振频率不同，在不同的发声方式下，特定频率的相对强弱也不同。于是我们只需要增强特定的频段就可以合成元音，这些较强的地方就被称为共振峰，即 ```Formant```，通常来说我们只关心前几个共振峰，因为高频通常会很弱，以下如果有需要会记 (fx) 为第 x 共振峰。

我实现程序其实有一个致命的 bug，因为它采用了直接合成特定强度的正弦波并叠加的方式，而不是先合成再用 FFT 增强特定频率的方式，我不知道这会不会有什么问题，因为我物理并不是很好。

话说回来，双元音的合成是平滑过渡共振峰的，尤其在日文中有所体现（均等过渡），不过这个项目不是合成日文的而是英文的，所以加入了可变参数控制两个音的比例。

为了让声音更加柔和，我加入了类似 ADSR 的发音强弱变化，这对于以下说明的辅音也是有用的。

### Consonants

辅音的发声原理就略显复杂了，大概分为好几类：

- 清辅音一类，只靠湍流，即高斯噪声，和各种音量变化，比如爆破音突然变化的音量，就可以合成。（这一类由于没有声带参与发声，只和特定部位的结构有关，所以没有共振峰，只有振动频率区间的区别）。
- 浊辅音一类，在清辅音的基础上加入一些类似共振峰的东西表示声带轻微振动即可。
- 鼻音边音一类，这个会有非常明显的共振峰，鼻音的共振峰非常明显，而边音的共振峰还会突然变化，这个令我很头疼，暂时没有想好怎么模拟。

对于音和音之间我还加入了一些交变参数使之变化不是很突然。

## Function introduction

编译要求在支持 C++17 以上的编译器进行（并不保证能在 **仅支持** C++14 的编译器中编译通过），编译命令在 Windows 下形如：

```g++ .\SoundTrans.cpp main.cpp -o main -O2 -Wall -std=c++17 -L"./fftw/" -lfftw3-3 -lfftw3f-3 -lfftw3l-3```

```main.cpp```，实现主要合成的功能，输入 BPM 和曲谱路径即可。

```SoundTrans.h/cpp``` 实现一些参数的转变，以及 FFT Filter，作为特定频率均衡器。

```Singer.h``` 规定 Singer 结构体的内容。

```Daisy.txt/high.txt/test.txt``` 给出了曲谱的格式示例。

```./Dictionary/``` 路径表示声库参数，具体如下，我将音符分为四类：

0. 清辅音，输入较强的部分
1. 元音，输入共振峰和对应强度
2. 浊辅音，输入同清辅音
3. 交变，双元音一类，输入交变的音，两者分别占用时长，以及中间交变段的占用百分比（占用各自的时长的百分比，这样可以保证短音被交变得更少）。

