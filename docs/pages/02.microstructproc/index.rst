.. _microstructproc:

.. title:: MPC

Microstructural profile covariance
============================================================

This module samples intracortical intensities from a microstructurally-sensitive contrast. This is achieved by constructing a series of equivolumetric surfaces between pial and white matter boundaries, yielding unique intracortical intensity profiles at each vertex of the native surface mesh. This approach has been previously applied in the `whole cortex <https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000284>`_, as well as in targeted structures such as the `insula <https://www.sciencedirect.com/science/article/pii/S1053811920303451>`_. By parcellating and cross-correlating nodal intensity profiles, this module profiles microstructural profile covariance (MPC) matrices describing similiarity in intracortical microstructure across the cortex.

Intracortical equivolumetric surfaces are generated using scripts from the `surface tools repository <https://github.com/kwagstyl/surface_tools>`_, available via GitHub.

.. figure:: sankey_mpc.png
   :align: center

-MPC
--------------------------------------------------------

.. admonition:: Prerequisites 🖐🏼

    You need to run ``-proc_structural``, ``-proc_freesurfer`` and ``-post_structural`` before this stage

.. figure:: MPC_01-14.gif
   :align: center

.. tabs::

    .. tab:: Processing steps

        - Compute boundary-based registration from microstructural imaging volume to FreeSurfer native space
        - Generate 16 equivolumetric surfaces between pial and white matter boundary previously defined from FreeSurfer. Surfaces closest to pial and white matter boundaries are then discarded to account for partial volume effects, resulting in 14 surfaces used for further analyses
        - Perform surface-based registration to fsaverage5 and conte69-32k templates
        - Average intensity profiles within parcels defined on the native FreeSurfer surface, while excluding outlier vertices
        - Perform partial correlation, controlling for parcel-wide mean profile, across all pairs of intensity profiles

    .. tab:: Usage

        **Terminal:**

        .. parsed-literal::
            $ mica-pipe **-sub** <subject_id> **-out** <outputDirectory> **-bids** <BIDS-directory> **-MPC** <options>

        **Docker command:**

        .. parsed-literal::
            $ docker -MPC

        **Optional arguments:**

        ``-MPC`` has several optional arguments:

        .. list-table::
            :widths: 75 750
            :header-rows: 1
            :class: tight-table

            * - **Optional argument**
              - **Description**
            * - ``-microstructural_img`` ``<path>``
              - Specifies an input image on which to sample intensities for the MPC analysis. You must specify this flag if your dataset does not include a qT1 image, or if your microstructurally-sensitive imaging contrast is not stored in the *rawdata* branch of the BIDS directory (for example, T1-weighted divided by T2-weighted derivative file). 
            * - ``-microstructural_lta`` ``<path>``
              - This option lets the user specify their own registration file to map the input image to native freesurfer space. The registration file must be in ``.lta`` format. If omitted, the registration will be performed in the script using `bbregister <https://surfer.nmr.mgh.harvard.edu/fswiki/bbregister/>`_.
            * - ``-microstructural_reg`` ``<path>``
              - Path to file which will be registered to native freesurfer space (e.g. ``./img_2reg.nii.gz``). This image can be different from the input provided to ``-microstructural_img``, but the two images must be in the same space!

    .. tab:: Outputs

        Directories created or populated by **-MPC**:

        .. parsed-literal::

            - <outputDirectory>/micapipe/<sub>/anat/surfaces/micro_profiles
            - <outputDirectory>/micapipe/<sub>/anat
            - <outputDirectory>/micapipe/<sub>/xfm

        Files generated by **-MPC**:

        .. parsed-literal::
            - Microstructural image in native FreeSurfer space:
               *<outputDirectory>/micapipe/<sub>/anat/<sub>_space-fsnative_desc-micro.nii.gz*

            - Boundary-based registration outputs from microstructural image to native FreeSurfer space:
               *<outputDirectory>/micapipe/<sub>/xfm/*
                   <sub>_from-micro_to-fsnative_bbr.dat
                   <sub>_from-micro_to-fsnative_bbr.dat.log
                   <sub>_from-micro_to-fsnative_bbr.dat.mincost
                   <sub>_from-micro_to-fsnative_bbr.dat.param
                   <sub>_from-micro_to-fsnative_bbr.dat.sum

            - Equivolumetric surface sampling ouputs stored in *<outputDirectory>/micapipe/anat/surfaces/micro_profiles*:

                - MPC json card:
                   *<sub>_MPC.json*

                - Native FreeSurfer space intensity at depth #:
                   *<sub>_space-fsnative_desc-<hemi>_MPC-<#>.mgh*

                - fsaverage5 template space intensity at depth #:
                   *<sub>_space-fsaverage5_desc-<hemi>_MPC-<#>.mgh*

                - conte69 template space intensity at depth #:
                   *<sub>_space-conte69-32k_desc-<hemi>_MPC-<#>.mgh*

                - Parcellated intensity profiles:
                   *<sub>_space-fsnative_atlas-<atlas>_desc-intensity_profiles.txt*

                - MPC matrices:
                   *<sub>_space-fsnative_atlas-<atlas>_desc-MPC.txt*
