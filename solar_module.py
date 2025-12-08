import pvlib

def calculateZenithAngles(coordinates, timedata, time_zone=None, UTC_offset=None):
    """ For given latitude, longitude, and set of datetimes, calculates the solar zenith angle. """
    latitude, longitude = coordinates
    
    #if time_zone is not None:
    #    timedata = timedata.dt.tz_localize(time_zone, ambiguous=False, nonexistent="shift_forward")
    #    values = pvlib.solarposition.get_solarposition(timedata, latitude, longitude)
    #    timedata = timedata.dt.tz_localize(None)
    #    return values
    if UTC_offset is not None:
        return pvlib.solarposition.get_solarposition(timedata - UTC_offset, latitude, longitude)
    else:
        return pvlib.solarposition.get_solarposition(timedata, latitude, longitude)
    
    
def calculateProjectedZenithAngles(coordinates, timedata, UTC_offset, zenith_offset, azimuth_offset):
    """ Should transpose the zenith angle to be relative to the plane of array
        WARNING: Does not appear to work correctly """
    solar_positions = calculateZenithAngles(coordinates, timedata, UTC_offset)
    return pvlib.shading.projected_solar_zenith_angle(solar_positions["zenith"], solar_positions["azimuth"], zenith_offset, azimuth_offset)




