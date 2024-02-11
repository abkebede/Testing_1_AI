import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.visualization import astropy_mpl_style
from astropy import units as u
from astroquery.simbad import Simbad

# Set the Astropy plotting style
plt.style.use(astropy_mpl_style)

# Define the function to plot the constellation outline and stars
def plot_constellation(constellation, query):
    # Query the Simbad star catalog for the constellation
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('ra(d)', 'dec(d)', 'flux(V)', 'sptype')
    stars = custom_simbad.query_region(constellation, radius=2 * u.deg)

    # Get the brightest stars in the constellation
    sorted_stars = sorted(stars, key=lambda star: star['FLUX_V'])
    brightest_stars = sorted_stars[:5]

    # Print the brightest stars and their magnitudes
    print('Brightest stars:')
    for star in brightest_stars:
        print(star['MAIN_ID'], star['FLUX_V'])

    # Get the coordinates of the constellation boundary
    boundary = SkyCoord.from_name(constellation).boundary()

    # Plot the boundary and stars
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlabel('Right Ascension (degrees)')
    ax.set_ylabel('Declination (degrees)')

    for segment in boundary:
        ra = segment.ra.deg
        dec = segment.dec.deg
        ax.plot(ra, dec, '-', color='grey')

    for star in stars:
        ra = star['RA_d']
        dec = star['DEC_d']
        mag = star['FLUX_V']
        size = np.sqrt(mag) / np.sqrt(4.83)
        ax.scatter(ra, dec, s=size * 100, edgecolors='none', color='white', alpha=0.5)

    # Set the title of the plot
    ax.set_title(f'Constellation {constellation} ({query})')

    # Show the plot
    plt.show()

# Example usage
constellation = 'Orion'
query = 'Betelgeuse'
plot_constellation(constellation, query)
