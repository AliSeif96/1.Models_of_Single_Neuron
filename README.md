# 1.Models_of_Single_Neuron
Simple Models of Single Neuron in Rodent Brains  
# WB
## Wang-Buzsaki Model of an Inhibitory Interneuron in Rat Hippocampus.

Wang and Buzsaki proposed a model of an inhibitory basket cell in rat hippocampus. Basket cells derive their name from the fact that the branches of their
axonal arbors form basket-like structures surrounding the cell bodies of othercells.

Two different classes of inhibitory basket cells are ubiquitous in the brain, the parvalbumin-positive (PV+) basket cells, which contain the protein parvalbumin, and the cholecystokinin-positive (CCK+) basket cells, which contain the hormone cholecystokinin. The PV+ basket cells are called fast-firing because they are capable of sustained high-frequency firing, and are known to play a central role in the generation of gamma frequency (30–80 Hz) oscillations. It is thought that gamma rhythms are important for sensory processing, attention, and working memory. The WB model is patterned after the fast-firing PV+ basket cells.

 The functions αx and βx are
 
<p align="center">
 <img src="https://github.com/aliseif321/1.Models_of_Single_Neuron/blob/main/1.____(WB)______Wang-Buzsaki%20Model%20Single%20Neuron/Book/2.png?raw=true" >
 </p>

Figure 5.3 in book shows a voltage trace with I = 0.75 μA/cm2.

<p align="center">
 <img src="https://github.com/aliseif321/1.Models_of_Single_Neuron/blob/main/1.____(WB)______Wang-Buzsaki%20Model%20Single%20Neuron/Book/1.png?raw=true" >
 </p>

The most striking difference between Figs. 5.3 is that the spike afterhyperpolarization, i.e., the hyperpolarization following an action potential, is far less deep in the WB model than in the RTM model. The difference between the lowest value of v and the firing threshold is about 15 mV in Fig. 5.3. This is in agreement with experimental results for fast-firing inhibitory interneurons; see references in. (See, however, also the voltage traces of cortical interneurons and pyramidal cells in Figs. 1C and 1E of. There the spike afterhyperpolarization is significantly more pronounced in the interneurons than in the pyramidal cells.) The spike afterhyperpolarization is less pronounced for the WB model than for the RTM model because the maximal conductance densities gNa and gK are smaller. Deeper afterhyperpolarization would be obtained if gNa and gK were raised (exercise 5), or h and n made slower. In fact, the Wang-Buzs´aki model as stated in included a scaling factor φ in front of the formulas for αh, βh, αn, and βn. Wang and Buzs´aki chose φ = 5. This choice is built into the equations as stated above. However, they pointed out that reducing φ, which amounts to reducing αh, βh, αn, and βn, i.e., to slowing down h and n, makes spike afterhyperpolarization more pronounced.


## C++
### Wang-Buzsaki model with Runge-Kutta 4th Order Method for one neuron




Topic: Wang-Buzsaki model with Runge-Kutta 4th Order Method for one neuron    Ali-Seif



Version Release 17.12 rev 11256



Date: 2/27/2021



Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler



MSI: PX60 6QD/ DDR4



Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM



<p align="center">
 <img src="https://github.com/aliseif321/1.Models_of_Single_Neuron/blob/main/1.____(WB)______Wang-Buzsaki%20Model%20Single%20Neuron/C++/Picture/pic1.png?raw=true" >
 </p>
