import os
from PIL import Image
import numpy as np
import anndata
import cairosvg
from colormap import rgb2hex
import squidpy as sq
import matplotlib.pyplot as plt

def svg_to_pixel_art(svg_path, output_png_dir, pixel_size=10):
    output_png_path = os.path.join(
       output_png_dir, "cairo.png")  # pixel_art1

    # Convert SVG to PNG
    cairosvg.svg2png(url=svg_path, write_to=output_png_path)

    # Open the image and resize to pixel art
    img = Image.open(output_png_path)
    width, height = img.size
    new_width = width // pixel_size
    new_height = height // pixel_size
    pixel_art = img.resize((new_width, new_height), Image.NEAREST)
    # Save the pixel art image
    pixel_art.save( os.path.join(
       output_png_dir, "pixel_art1.png") )

    return pixel_art

def pixel_art_to_anndata(pixel_art):
    # Convert pixel art to numpy array
    pixel_array = np.array(pixel_art)

    # Reshape to (n_pixels, 3) for RGB
    n_pixels = pixel_array.shape[0] * pixel_array.shape[1]
    colors = pixel_array.reshape(n_pixels, -1)

    # Create x, y coordinates for each pixel
    height, width = pixel_array.shape[:2]
    x_coords, y_coords = np.meshgrid(np.arange(width), np.arange(height))
    coords = np.vstack((x_coords.ravel(), y_coords.ravel())).T

    # Create anndata object
    adata = anndata.AnnData(X=colors)
    # Assign coordinates to obsm
    adata.obsm['spatial'] = coords

    # and the spatial params for plotting purposes
    adata.uns['spatial'] = {
        'MY_SLICE7_NAME': {'scalefactors':
                               {'fiducial_diameter_fullres': 100,
                                'spot_diameter_fullres': 1,
                                'tissue_hires_scalef': 0.5,
                                'tissue_lowres_scalef': 0.5}}}

    return adata


def recover_colors(adata):
    tmp = adata.X.T   #.to_dense()
    colors_hex = [rgb2hex(k, i, j) for (k, i, j) in zip(tmp[0], tmp[1], tmp[2])]

    adata.obs['group_truth'] = colors_hex

    return adata



def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')
    print(rgb2hex(255,255,255, normalised=False))
    # Example usage
    svg_path = os.path.join(
        "svg_fig", "image1.svg"
    )
    output_png_dir = "png_out"

    pixel_art = svg_to_pixel_art(svg_path, output_png_dir, pixel_size=10)
    adata = pixel_art_to_anndata(pixel_art)
    adata = recover_colors(adata)

    # Save the anndata object
    adata.write("pixel_art.h5ad")
    print(adata.to_df())
    sq.pl.spatial_scatter(adata, color='group_truth',
    palette=adata.obs['group_truth'].tolist(),
    img=False)
    plt.show()




# See PyCharm help at https://www.jetbrains.com/help/pycharm/
