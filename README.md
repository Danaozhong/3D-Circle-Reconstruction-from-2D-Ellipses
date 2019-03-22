# 3D Circle Reconstruction from 2D Ellipses on a Camera Image

If you have serveral camera pictures of a 3D circle taken from different orientations, there is a closed-form method to reconstruct the orientation and diameter of the 3D circle with respect to the different camera positions.

Basically you need to know you camera parameters (intrinsic parameters), and the relative positions of your camera, and the ellipsis 2D data on the camera image (can be retrieved using a canny edge detector, followed by RANSAC / least square fitting of ellipses).

The result is the position, orientation and diameter of the 3D circle.

The implementation is based on the research paper below:

Soheilian, B.; Brédif, M., “Multi-view 3D Circular Target Reconstruction with Uncertainty Analysis”, ISPRS Annals of the Photogrammetry, Remote Sensing and Spatial Information Sciences, Vol. 2, Issue 3, pp 143-148, 2014


All technical details described here:
http://texelography.com/2019/01/07/3d-circle-reconstruction-from-ellipses/
