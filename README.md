# Grid Class

Cartesian Cell-based Anisotropic Refinement (CCAR) Grid.  Performs adaptive anisotropic mesh refinement for flow simulation.
Visualizations included the grid itself, with or without cell centers marked, and the velocity field.  

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
@article{berger1989local,
	 Author = {Berger, Marsha J and Colella, Phillip},
	 Annote = {},
	 Journal = {Journal of Computational Physics},
	 Keywords = {gridding, adaptive mesh refinement},
	 Number = {1},
	 Pages = {64--84},
	 Publisher = {Elsevier},
	 Title = {Local adaptive mesh refinement for shock hydrodynamics},
	 Volume = {82},
	 Year = {1989}
}

@article{berger1984adaptive,
	 Author = {Berger, Marsha J and Oliger, Joseph},
	 Annote = {},
	 Keywords = {gridding, adaptive mesh refinement},
	 Journal = {Journal of Computational Physics},
	 Number = {3},
	 Pages = {484--512},
	 Publisher = {Elsevier},
	 Title = {Adaptive mesh refinement for hyperbolic partial differential equations},
	 Volume = {53},
	 Year = {1984}
}

@incollection{davis2017adaptive,
  author={Davis, P},
  title={Adaptive Mesh Refinement: An Essential Ingredient in Computational Science.},
  booktitle={SIAM News},
  year={2017}, 
  Keywords = {AMR, Gridding},
  Annote={}
}

@inproceedings{nilsson2005novel,
  author = {Nilsson, J. and Gerritsen, M. and Younis, R. and others},
  title = {A novel adaptive anisotropic grid framework for efficient reservoir simulation},
  booktitle = {SPE reservoir simulation symposium},
  organization = {Society of Petroleum Engineers},
  year = {2005},
  keywords = {CCAR, grid adaptivity, anisotropic},
  annote = {Development of CCAR for anisotropy.}
}
``` 
