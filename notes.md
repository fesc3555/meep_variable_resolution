This gives a little overview of the procedure used to gain variable resolution using <span style="font-variant:small-caps;">meep</span> as a library in <span style="font-variant:small-caps;">C++</span>. This document was intended purely as a note for myself and should be read as such. I programmed a couple of examples that should be grouped together in archives with this document. The code is probably neither fast nor well written, but it served its purpose for me.

The geometry transformation is sketched in figure \[fig:sketch\].

<embed src="skizze.eps" style="width:50.0%" />

The coordinate stretching can be absorbed into a transformation of *ε*, *μ*, *ρ*, **J** and the fields in Maxwells Equations see for example S. Johnsons’s notes on a MIT course, “Coordinate Transformation & Invariance in Electromagnetism”.

Implementation Details
======================

Define a grid volume as usual, using the maximal resolution of the central region. Calculate the appropriate scaled size beforehand and reduce the size of the outer regions accordingly when defining the volume.

Define a class that handles the different resolution domains, including functions that give the scale factors for the *χ*<sub>1</sub> + 1, i.e. *ε* or *μ* respectively, transformation and the **J** transformation. Note that *χ*<sub>1</sub> + 1 now definitely is a (diagonal) tensor in the domains, where the stretched geometry has rectangular pixels, so the stretching factors depend on the field component. The input for these member functions is therefore a position vector and a component direction.

Next, define a material property class. Because we want to specify tensor-valued *χ*<sub>1</sub> + 1, we directly have to override the virtual `eff_chi1_inv` function. This function also handles anisotropic averaging, so this also has to be implemented. **You might be tempted to ignore the anisotropic averaging, but this is very unwise, since the resolution-domain boundaries will cause unphysical reflections if you forgo the averaging.** Cascading the resolution, i.e. using more domains in order to reduce the difference between adjacent domains, even worsens reflection in my examples (I just made one large jump). Tests show, however, that you can significantly rise the tolerances and lower the maximum iterations of this (originally fall-back) procedure without compromising the physics very much. So, for now, just copy and paste the anisotropic avergaing code into our function, together with the `normal_vector` (renamed to `normal_vector_anis` in the exmaple) function and the `sphere_point` helper function (import sphere-quad.h for it to work). In `eff_chi1_inv`, now replace the line where it calls `chi1p1` with the stretching function times `chi1p1`. Do the same in the `normal_vector_anis` routine, adding an additional argument for the component direction in the function header and changing the function call in `eff_chi1_inv` accordingly. I actually have not thought this through, but this way of calculating a gradient seems weird, multiplying a gradient by a tensor component. It works, however, for now. Additionally, reimplement `has_mu` and let it answer `true` because of the transformation.

Now for the source transformation. In general, to get plane waves, a volume source is necessary (“volume” in this case being a hypersurface). If the volume extends over several resolution domains, the source amplitude in the different domains has to be adjusted. While this is possible in meep by giving `add_volume_source` a function pointer, it is not convenient for the given implementation, i.e. using a class that handles the resolution domains. We want to give `add_volume_source` a member function of our variable-resolution handler class instead of a function pointer, two totally different things in <span style="font-variant:small-caps;">C++</span>. Converting the latter into the former is nearly impossible and, if done, unreadably ugly. So instead, I changed the meep library source files and recompiled. Specifically, I changed the function header of `add_volume_source` so it takes a `std::function<complex<double>(const vec &)>` object instead of a function pointer. This assures backwards compatibilty (the function object can be initialised with a function pointer) and every member function can be cast into this std::function with std::bind, at least if the compiler is not at least a decade outdated and knows its <span style="font-variant:small-caps;">C++11</span>. The compiler reminded me to make the same replacement in the struct `source_vol_chunkloop_data`. Important at this point is that the source amplitudes are, for whatever reason, given in a grid volume centered coordinate system. Even if normally the lower left corner is at (0, 0, 0), now that vector points to the middle of the system, apparently. Correcting for that is, however, easy.

Scattering Spectra
==================

Scattering spectra are usually calculated with `dft_flux_box`-es. This method calculates the frequency-dependent energy flux through each point of a surface enclosing the scatterer. First, a calculation without the scatterer is done to get the flux caused by the incoming wave. Afterwards, another calculation including the scatterer has to done and the flux of the incident wave from the former calculation has to be substracted at each point. Finally, the energy flux is integrated over the surface, giving the overall scattering intensity, which afterwards has to be normalized with regard to the frequency-dependent excitation intensity, i.e. the Fourier-transform of the incoming pulse.

The question in the context of a variable-resolution geometry is, whether the flux through a surface is invariant under the given coordinate transform. We can show that it is indeed, meaning that the routine described above can be applied in the variable-resolution calculations without any modifications.

We want to show the following equality:
$$\\begin{aligned}
\\iint\_\\Omega \\textbf{S}(\\omega) \\cdot \\mathrm{d}\\textbf{A}\_\\Omega =
\\iint\_{\\Omega'} \\textbf{S}'(\\omega) \\cdot \\mathrm{d}\\textbf{A}'\_{\\Omega'} \\quad,
\\label{eq:toshow}\\end{aligned}$$
 with the surface *Ω*, the time-averaged Poynting vector $\\textbf{S}(\\omega) = \\frac{1}{2} \\textbf{E}(\\omega) \\times \\textbf{H}^\*(\\omega)$, the vector area element (normal vector times scalar area element) on the the surface *d***A**<sub>*Ω*</sub>. Primed quantities are the transformed quantities according to S. Johnson’s reference. The surface *Ω* is given by a parametrization **φ**(*u*, *v*). The vector area element is
$$\\begin{aligned}
\\mathrm{d}\\textbf{A}\_\\Omega &= \\boldsymbol{\\varphi}\_u \\times \\boldsymbol{\\varphi}\_v  \\, \\mathrm{d}u \\, \\mathrm{d}v \\\\
\\boldsymbol{\\varphi}\_u &= \\frac{\\partial \\boldsymbol{\\varphi}}{\\partial u} \\nonumber\\end{aligned}$$
 Changing to index notation and employing Einstein’s sum convention, we get
$$\\begin{aligned}
2\\iint\_\\Omega \\textbf{S}(\\omega) \\cdot \\mathrm{d}\\textbf{A}\_\\Omega &=
\\iint\_\\Omega \\varepsilon\_{abc} E\_a H\_b^\* \\mathrm{d}A\_c\\nonumber  \\\\
&=\\iint\_\\Omega \\varepsilon\_{abc} E\_a H\_b^\* \\,\\varepsilon\_{mnc} {\\varphi\_m}\_u {\\varphi\_n}\_v  \\, \\mathrm{d}u \\, \\mathrm{d}v \\nonumber \\\\
&=\\iint\_{\\Omega'} \\varepsilon\_{abc} \\mathcal{J}\_{i a} E'\_i \\mathcal{J}\_{j b} H'^\*\_j \\,
\\varepsilon\_{mnc} \\mathcal{J}^{-1}\_{mx} {\\varphi'\_x}\_u \\mathcal{J}^{-1}\_{ny} {\\varphi'\_y}\_v  \\, \\mathrm{d}u' \\, \\mathrm{d}v'\\nonumber \\\\
&= \\iint\_{\\Omega'} \\varepsilon\_{abc} \\mathcal{J}\_{i a} E'\_i \\mathcal{J}\_{j b} H'^\*\_j \\,\\, \\delta\_{cd} \\,\\, 
\\varepsilon\_{mnd} \\mathcal{J}^{-1}\_{mx} {\\varphi'\_x}\_u \\mathcal{J}^{-1}\_{ny} {\\varphi'\_y}\_v  \\, \\mathrm{d}u' \\, \\mathrm{d}v'\\nonumber \\\\
&= \\iint\_{\\Omega'} \\varepsilon\_{abc} \\mathcal{J}\_{i a} E'\_i \\mathcal{J}\_{j b} H'^\*\_j \\,\\,
\\delta\_{dc} \\,\\, 
\\varepsilon\_{mnd} \\mathcal{J}^{-1}\_{mx} {\\varphi'\_x}\_u \\mathcal{J}^{-1}\_{ny} {\\varphi'\_y}\_v  \\, \\mathrm{d}u' \\, \\mathrm{d}v'\\nonumber \\\\
&= \\iint\_{\\Omega'} \\varepsilon\_{abc} \\mathcal{J}\_{i a} E'\_i \\mathcal{J}\_{j b} H'^\*\_j \\,\\,
\\mathcal{J}^{-1}\_{d \\alpha} \\mathcal{J}\_{\\alpha c} \\,\\, 
\\varepsilon\_{mnd} \\mathcal{J}^{-1}\_{mx} {\\varphi'\_x}\_u \\mathcal{J}^{-1}\_{ny} {\\varphi'\_y}\_v  \\, \\mathrm{d}u' \\, \\mathrm{d}v'\\nonumber \\\\
&= \\iint\_{\\Omega'} \\varepsilon\_{abc} \\mathcal{J}\_{i a} \\mathcal{J}\_{j b} \\mathcal{J}\_{\\alpha c} E'\_i H'^\*\_j \\,\\, 
\\varepsilon\_{mnd} \\mathcal{J}^{-1}\_{mx} \\mathcal{J}^{-1}\_{ny} \\mathcal{J}^{-1}\_{d \\alpha}{\\varphi'\_x}\_u {\\varphi'\_y}\_v  \\, \\mathrm{d}u' \\, \\mathrm{d}v'\\nonumber \\\\
&= \\iint\_{\\Omega'} \\varepsilon\_{ij\\alpha} \\mathrm{det}(\\mathcal{J}) E'\_i H'^\*\_j \\,\\, 
\\varepsilon\_{xy\\alpha} \\mathrm{det}(\\mathcal{J}^{-1}){\\varphi'\_x}\_u {\\varphi'\_y}\_v  \\, \\mathrm{d}u' \\, \\mathrm{d}v' \\nonumber \\\\
&= \\iint\_{\\Omega'} \\varepsilon\_{ij\\alpha} E'\_i H'^\*\_j \\,\\, \\varepsilon\_{xy\\alpha} {\\varphi'\_x}\_u {\\varphi'\_y}\_v  \\, \\mathrm{d}u' \\, \\mathrm{d}v' \\quad , \\end{aligned}$$
 which is just Eq. \[eq:toshow\]. Here, *δ*<sub>*c**d*</sub> denotes the Kronecker-delta. In the process, we used the fact that
$$\\begin{aligned}
{\\varphi\_m}\_u = \\frac{\\partial \\varphi\_m}{\\partial u} = \\frac{\\partial \\varphi\_m}{\\partial \\varphi'\_x} \\frac{\\partial \\varphi'\_x}{\\partial u} = \\mathcal{J}^{-1}\_{mx} {\\varphi'\_x}\_u\\end{aligned}$$
 as well as *ε*<sub>*a**b**c*</sub>𝒥<sub>*i**a*</sub>𝒥<sub>*j**b*</sub>𝒥<sub>*k**c*</sub> = *ε*<sub>*i**j**k*</sub>*d**e**t*(𝒥), which is implied by the definition of the determinant of 3 × 3 matrices *ε*<sub>*a**b**c*</sub>𝒥<sub>1*a*</sub>𝒥<sub>2*b*</sub>𝒥<sub>3*c*</sub> = *d**e**t*(𝒥), and the well known relations *d**e**t*(*A*)=*d**e**t*(*A*<sup>*T*</sup>)=1/*d**e**t*(*A*<sup>−1</sup>).

Example files in this folder
============================

The example files contain:

-   a 2D and a 3D example geometry with a dielectric sphere (cylinder) in the middle, pmls and point source.

-   a 2D example with a dielectric cylinder illuminated by a plane wave source. The folder also contains the changed meep source files.

-   a 2D example that calculates the scattering of a gold cylinder using `dft_flux_box`.

I did not include any evaluation routines.
