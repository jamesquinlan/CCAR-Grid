# Grid Class

Cartesian Cell-based Anisotropic Refinement (CCAR) Grid.  Performs adaptive anisotropic refinement for accurate and robust flow simulation.

```
MATLAB Version: 9.10.0.1602886 (R2021a)
Quinlan, James
```
 


## Gridi

Gridi, (Grid with an "i"), is an attempt to improve the original grid class.  Improvements include:

1. Removing edge and corner cells ("fake cells"). Perhaps keep the grid skeleton. 
2. Create function to get cell ID (and FaceID).
3. Improve kron creation.  Look for idea, placed in comments, about generating without kron. 
4. Consider private properties and functions in some cases.  
5. Consider creating @classFolder.
6. Consider creating small independent functions for constants, etc.  


## References

```
@incollection{davis2017adaptive,
  title={Adaptive Mesh Refinement: An Essential Ingredient in Computational Science.},
  author={Davis, P},
  booktitle={SIAM News},
  year={2017}, 
  Keywords = {AMR, Gridding},
  Annote={}
}

@inproceedings{nilsson2005novel,
	Author = {Nilsson, J. and Gerritsen, M. and Younis, R. and others},
	Keywords = {CCAR, grid adaptivity, anisotropic},
	Annote = {Development of CCAR for anisotropy. },
	Booktitle = {SPE reservoir simulation symposium},
	Organization = {Society of Petroleum Engineers},
	Title = {A novel adaptive anisotropic grid framework for efficient reservoir simulation},
	Year = {2005}}
``` 
