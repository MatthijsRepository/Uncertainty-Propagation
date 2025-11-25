import pvlib

def calculateZenithAngles(coordinates, timedata, UTC_offset):
    latitude, longitude = coordinates
    
    if UTC_offset is not None:
        return pvlib.solarposition.get_solarposition(timedata - UTC_offset, latitude, longitude)
    else:
        return pvlib.solarposition.get_solarposition(timedata, latitude, longitude)
    
    
def calculateProjectedZenithAngles(coordinates, timedata, UTC_offset, zenith_offset, azimuth_offset):
    solar_positions = calculateZenithAngles(coordinates, timedata, UTC_offset)
    return pvlib.shading.projected_solar_zenith_angle(solar_positions["zenith"], solar_positions["azimuth"], zenith_offset, azimuth_offset)




