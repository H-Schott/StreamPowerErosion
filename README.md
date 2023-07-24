# Stream Power Erosion source code

This is the release code for the paper *Large-scale terrain authoring through interactive erosion simulation*, to be published in Transactions on Graphics.

A preprint is available on the french platform HAL : [https://hal.science/hal-04049125v1](https://hal.science/hal-04049125v1)

The user paints an **uplift** map, which specifies the local elevation speed of the ground. Fluvial erosion is continuously simulated through the **Stream Power Equation** from geomorphology. The equilibrium between those two forces results in a dendritic mountainous terrain with a correct outflowing water network. Mountainous ridges appear where high uplift values have been placed, which allows to interactively author the terrain.

To cite this code or the related article :
```tex
@article{Schott2023,
author = {Schott, Hugo and Paris, Axel and Fournier, Lucie and Gu\'{e}rin, Eric and Galin, Eric},
title = {Large-Scale Terrain Authoring through Interactive Erosion Simulation},
year = {2023},
journal = {ACM Trans. Graph.}
}
```


## Compile

Clone the repository and its submodules using:
```
git clone --recursive https://github.com/H-Schott/StreamPowerErosion.git
```

Can be compile on Windows or Linux, using the CMake file.


## Control

**ctrl + mouse wheel** to zoom in and out.  
**ctrl + left click** to add uplift and place mountains.