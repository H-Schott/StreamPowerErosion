# Stream Power Erosion source code

This is the release code for the paper *Large-scale terrain authoring through interactive erosion simulation*, to be published in Transactions on Graphics.

A preprint is available on the french platform HAL : [https://hal.science/hal-04049125v1](https://hal.science/hal-04049125v1)

The user paints an **uplift** map, which specifies the local elevation speed of the ground. Fluvial erosion is continuously simulated through the **Stream Power Equation** from geomorphology. The equilibrium between those two forces results in a dendritic mountainous terrain with a correct outflowing water network. Mountainous ridges appear where high uplift values have been placed, which allows to interactively author the terrain.

To cite this code or the related article :
```tex
@unpublished{schott:hal-04049125,
  TITLE = {{Large-scale terrain authoring through interactive erosion simulation}},
  AUTHOR = {Schott, Hugo and Paris, Axel and Fournier, Lucie and Gu{\'e}rin, Eric and Galin, Eric},
  URL = {https://hal.science/hal-04049125},
  NOTE = {working paper or preprint},
  YEAR = {2023},
  MONTH = Mar,
  KEYWORDS = {Erosion simulation ; Landscapes},
  PDF = {https://hal.science/hal-04049125/file/2022-uplift-author.pdf},
  HAL_ID = {hal-04049125},
  HAL_VERSION = {v1},
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