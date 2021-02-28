# Erisir
## Erisir Model of an Inhibitory Interneuron in Mouse Cortex.

Erisir et al. [45] proposed a model of an inhibitory interneuron in mouse somatosensory cortex. With minor modifications discussed in detail in [17], the model takes the same form as the RTM and WB models, except that the potassium conductance is gKn2, not gKn4; the significance of this difference will be discussed at the end of this section. For the constants. The functions αx and βx are

 The functions αx and βx are
 
<p align="center">
 <img src="https://github.com/aliseif321/1.Models_of_Single_Neuron/blob/main/3.____(Erisir)______ErisirModel%20Single%20Neuron/Book/Untitled.png?raw=true" >
 </p>

Figure 5.4 shows a voltage trace with I = 7 μA/cm2.

<p align="center">
 <img src="https://github.com/aliseif321/1.Models_of_Single_Neuron/blob/main/3.____(Erisir)______ErisirModel%20Single%20Neuron/Book/Untitled2.png?raw=true" >
 </p>

Note that gNa and gK are quite large in the Erisir model, even larger than in the RTM model. As a result, the voltage rises almost to vNa during an action potential, and falls almost to vK immediately following an action potential. The leak conductance density gL is large as well.
What is the significance of taking the potassium conductance to be gKn2, not gKn4? The main answer is that it does not appear to matter very much; compare Fig. 5.5, where we have used gKn4 instead of gKn2, with Fig. 5.4. In detail, using gKn2 instead of gKn4 has the following effects, which one can see when comparing Figs. 5.4 and 5.5.

1. As n rises to values near 1 during a spike, the potassium conductance responds more rapidly when the exponent is 2, not 4. Therefore the spike termination mechanism becomes faster, and the spikes become narrower.
2. As n falls to values near 0 following a spike, the potassium conductance follows less rapidly when the exponent is 2, not 4. This has the effect that the hyperpolarization following a spike is deeper.
3. Surprisingly, even though the potassium current is hyperpolarizing, and gKn2 is greater than gKn4 for 0 <n< 1, firing is slightly faster with gKn2 than with gKn4. In essence, this is explained by the fact that the narrower action potentials in Fig. 5.4 allow less time for deep inactivation of the sodium current.

<p align="center">
 <img src="https://github.com/aliseif321/1.Models_of_Single_Neuron/blob/main/3.____(Erisir)______ErisirModel%20Single%20Neuron/Book/Untitled3.png?raw=true" >
 </p>

## C++
### Erisir model with midle point Method for one neuron




Topic: Erisir model with midle point Method for one neuron    Ali-Seif



Version Release 17.12 rev 11256



Date: 3/2/2021



Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler



MSI: PX60 6QD/ DDR4



Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM



<p align="center">
 <img src="https://github.com/aliseif321/1.Models_of_Single_Neuron/blob/main/3.____(Erisir)______ErisirModel%20Single%20Neuron/C++/Picture/Untitled.png?raw=true" >
 </p>

<p align="center">
 <img src="https://github.com/aliseif321/1.Models_of_Single_Neuron/blob/main/3.____(Erisir)______ErisirModel%20Single%20Neuron/C++/Picture/Untitled4.png?raw=true" >
 </p>
