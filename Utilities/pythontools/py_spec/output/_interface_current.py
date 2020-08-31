########################################
# interface_current.py: functions that calculate the interface current
#

def interface_current(self, interface_number):
    """
    returns the current on the interface given by interface_number.
    Not working yet

    We are calculating a loop integral over the difference in magnetic field on the surface.
    All higher fourier harmonics fall out, so we only need to take the difference between the
    zeroth components just outside and in.
      Arguments:
        self: An instance of the py_self.self class run on the output.h5 file.
        interface_number: the number of the interface on which the current is to be calculated.
                          STARTS AT 1!
    """
    import numpy as np
    Ntor = self.physics.Ntor
    # fall out.
    outer_volume = self.output.Btemn[0:Ntor, 0, interface_number]
    inner_volume = self.output.Btemn[0:Ntor, 1, interface_number-1]

    current = np.sum(outer_volume-inner_volume) * np.pi * 2  # check if idl's pi2 is pi squared or just two pi

    return current


def interface_current_total(self):
    """
    NOT TESTED!
    returns the current on the interface given by interface_number.
    Not working yet

    We are calculating a loop integral over the difference in magnetic field on the surface.
    All higher fourier harmonics fall out, so we only need to take the difference between the
    zeroth components just outside and in.
      Arguments:
        self: An instance of the py_self.self class run on the output.h5 file.
        interface_number: the number of the interface on which the current is to be calculated.
    """
    import numpy as np
    number_of_interfaces = self.physics.Nvol - 1 # There is one less interface than there are volumes
    current = 0.

    for interface_number in range(1,number_of_interfaces):  #Not doing the last interface (range stops 1 early)
        current += interface_current(self, interface_number)

    return current





#for ivol = 0, Ndof-1 do values[ivol] = [ total( lBtemn[0:Ntor,0,ivol+1]-lBtemn[0:Ntor,1,ivol] ) ] * pi2
