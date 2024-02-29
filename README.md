# Coarse-Frequency-Translation
Coarse frequency translation of incoming signal with using FFT and required preprocessing steps. 

<br> This project is developed in the scope of Software Defined Communications lecture. The aim of the project is to estimate the coarse carrier frequency by using FFT and required preprocessing steps. The result of this project is preprocessing step of Phase-Locked-Loop. The functions of this project  is obtained from C. R. Johnson Jr., W. A. Sethares, and A. G. Klein. Software receiver design: build your own digital communication system in five easy steps. Cambridge University Press, 2011

<br>Contiributors:
<br> Yahya Ekin
<br> Özün Yiğit Bayram

<br>We have two version of our code.
<br>1-CoarseFreqTrans_SingleTime.m
<br>2-CoarseFreqTrans_MultipleTimes.m

<br>1-CoarseFreqTrans_SingleTime.m: Calculates Normalized Symbol Error Rate for one time for given freuency offset. Plots all the figures on everye step that given in presentation slide.

<br>2-CoarseFreqTrans_MultipleTimes.m: This version plots Normalized Symbol Error rate for given band of frequencies. Also, calculates error per offset for trial times. All the frequency offset band and trial time can be adjusted.

<br>Since the code is based on functions of the book, it works in matlab files folder provided by the book.
