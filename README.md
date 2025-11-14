# Tricritical Ising Conformal Defects

We will calculate some of the conformal defects that the tricritical Ising CFT admits in order to identify some stationary points under RG flow. This is a place to store the code needed for the relevant calculations. 



## How do I run this?

You need `sage`. To get it [check this out](https://doc.sagemath.org/html/en/installation/source.html). 

Once you have sage installed, we have build a python module with tools for CFT that you can use! In the future this would be instaleld via pip by
```
sage -python -m pip install cftpy
```
but for now, if you want to run the notebooks and edit everything, a local version of the module is in this repository under `./cftpy`. You can install it locally via
```
save -python -m pip install -e ./cftpy
```
the `-e` flag stands for editable, so if you `git pull` our updates sage will automatically use the updated module.


