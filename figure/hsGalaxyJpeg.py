"""@file hs_galaxy_image.py
Utilities to deal with galaxy image in JPEG or PNG format.
Including functions to scale, rotate, transform, and filter galaxy images to
study their morphology
"""

import numpy as np
import subprocess as sub
from os import path

def find_image_center(xsize,ysize):
    """
    Find the pixel coordinates of the image center
    """

    if xsize % 2 == 0:
        cen_x = (xsize+1)/2.0
    else:
        cen_x = xsize/2.0

    if ysize % 2 == 0:
        cen_y = (ysize+1)/2.0
    else:
        cen_y = ysize/2.0

    return (cen_x, cen_y)


def file_name_attach(input_file, suffix, dir_out=None):
    """
    Attach a suffix to the name of a file
    """

    # Split the path of the file
    dir_now = path.dirname(input_file)
    base_now = path.basename(input_file)
    pre, ext = path.splitext(base_now)
    pre = pre.strip()
    ext = ext.strip()

    # Suffix to attach
    suf = suffix.strip()

    # New file name
    base_new = pre + "_" + suf + ext

    # Location and name of the output image
    if dir_out is None:
        output_file = path.join(dir_now, base_new)
    else:
        output_file = path.join(dir_out.strip(), base_new)

    return output_file


def magick_imagesize(file_name):
    """
    Use the identify@ImageMagick to get the dimensions of the images
    in unit of pixel
    """
    proc = sub.Popen('identify ' + file_name, stdout=sub.PIPE, shell=True)
    size = proc.stdout.read().split()[2].split("x")

    xsize, ysize = int(size[0]), int(size[1])

    return (xsize, ysize)


def is_pixel_inside(dim, coord):
    """
    Check if the pixel coordinate is inside the image
    """

    if (len(dim)<2) or (len(coord)<2):
        raise Exception("Dimensions should be >= 2! Check!")
    if (0<=coord[0]<=dim[0]) and (0<=coord[1]<=dim[1]):
        return True
    else:
        return False


def magick_rotate(input_file, ang_rot, suffix=None, center=None,
                  bg_color=None, dir_out=None):
    """
    Use the SRT distortion function in convert@ImageMagick to rotate an image by
    certain angle, and save the output
    """

    # Check the existence of the file
    if not path.exists(input_file):
        raise Exception("Can not find the input image:" + input_file + '!!')

    # Color for the virtual pixels
    if bg_color is None:
        bg_color = 'black'
    else:
        bg_color = str(bg_color)

    # Define the name of the output file
    ang_str = "{:.0f}".format(ang_rot)
    if suffix is None:
        suffix = ang_str.strip()
    else:
        suffix = suffix.strip()

    # Get the output file name
    output_file = file_name_attach(input_file, suffix, dir_out=dir_out)

    # Check the center for the rotation
    xsize, ysize=magick_imagesize(input_file)
    if center is None:
        cen_x, cen_y = find_image_center(xsize,ysize)
    else:
        cen_x, cen_y = center[0], center[1]
        if not is_pixel_inside([xsize,ysize],[cen_x,cen_y]):
            raise Exception("The rotation center should be inside the image!!")
    cen_x, cen_y = str(cen_x).strip(), str(cen_y).strip()

    # Call ImageMagick to do the rotation
    convert = '/opt/local/bin/convert'
    cmd = convert + ' ' + input_file + ' -virtual-pixel ' + bg_color + \
          ' -interpolate bicubic -distort SRT "' + cen_x + ',' + cen_y + \
            ',1.0,1.0,' + ang_str.strip() + ',' + cen_x + ',' + cen_y + '" ' + \
           output_file
    proc=sub.call(cmd, shell=True)
    if proc is not 0:
        raise Exception("Rotation is not done properly! Check!")

    return output_file

