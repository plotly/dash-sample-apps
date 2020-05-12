import aerosandbox as asb
from aerosandbox.library.airfoils import e216
import numpy as np
import casadi as cas
import copy

naca0008 = asb.Airfoil("naca0008")


def make_airplane(
    n_booms, wing_span,
):
    # n_booms = 3

    # wing
    # wing_span = 37.126
    wing_root_chord = 2.316
    wing_x_quarter_chord = -0.1

    # hstab
    hstab_span = 2.867
    hstab_chord = 1.085
    hstab_twist_angle = -7

    # vstab
    vstab_span = 2.397
    vstab_chord = 1.134

    # fuselage
    boom_length = 6.181
    nose_length = 1.5
    fuse_diameter = 0.6
    boom_diameter = 0.2

    wing = asb.Wing(
        name="Main Wing",
        # x_le=-0.05 * wing_root_chord,  # Coordinates of the wing's leading edge # TODO make this a free parameter?
        x_le=wing_x_quarter_chord,  # Coordinates of the wing's leading edge # TODO make this a free parameter?
        y_le=0,  # Coordinates of the wing's leading edge
        z_le=0,  # Coordinates of the wing's leading edge
        symmetric=True,
        xsecs=[  # The wing's cross ("X") sections
            asb.WingXSec(  # Root
                x_le=-wing_root_chord / 4,
                # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                y_le=0,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                z_le=0,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                chord=wing_root_chord,
                twist=0,  # degrees
                airfoil=e216,  # Airfoils are blended between a given XSec and the next one.
                control_surface_type="symmetric",
                # Flap # Control surfaces are applied between a given XSec and the next one.
                control_surface_deflection=0,  # degrees
                spanwise_panels=30,
            ),
            asb.WingXSec(  # Tip
                x_le=-wing_root_chord * 0.5 / 4,
                y_le=wing_span / 2,
                z_le=0,  # wing_span / 2 * cas.pi / 180 * 5,
                chord=wing_root_chord * 0.5,
                twist=0,
                airfoil=e216,
            ),
        ],
    )
    hstab = asb.Wing(
        name="Horizontal Stabilizer",
        x_le=boom_length
        - vstab_chord * 0.75
        - hstab_chord,  # Coordinates of the wing's leading edge
        y_le=0,  # Coordinates of the wing's leading edge
        z_le=0.1,  # Coordinates of the wing's leading edge
        symmetric=True,
        xsecs=[  # The wing's cross ("X") sections
            asb.WingXSec(  # Root
                x_le=0,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                y_le=0,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                z_le=0,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                chord=hstab_chord,
                twist=-3,  # degrees # TODO fix
                airfoil=naca0008,  # Airfoils are blended between a given XSec and the next one.
                control_surface_type="symmetric",
                # Flap # Control surfaces are applied between a given XSec and the next one.
                control_surface_deflection=0,  # degrees
                spanwise_panels=8,
            ),
            asb.WingXSec(  # Tip
                x_le=0,
                y_le=hstab_span / 2,
                z_le=0,
                chord=hstab_chord,
                twist=-3,  # TODO fix
                airfoil=naca0008,
            ),
        ],
    )
    vstab = asb.Wing(
        name="Vertical Stabilizer",
        x_le=boom_length - vstab_chord * 0.75,  # Coordinates of the wing's leading edge
        y_le=0,  # Coordinates of the wing's leading edge
        z_le=-vstab_span / 2
        + vstab_span * 0.15,  # Coordinates of the wing's leading edge
        symmetric=False,
        xsecs=[  # The wing's cross ("X") sections
            asb.WingXSec(  # Root
                x_le=0,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                y_le=0,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                z_le=0,  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                chord=vstab_chord,
                twist=0,  # degrees
                airfoil=naca0008,  # Airfoils are blended between a given XSec and the next one.
                control_surface_type="symmetric",
                # Flap # Control surfaces are applied between a given XSec and the next one.
                control_surface_deflection=0,  # degrees
                spanwise_panels=8,
            ),
            asb.WingXSec(  # Tip
                x_le=0,
                y_le=0,
                z_le=vstab_span,
                chord=vstab_chord,
                twist=0,
                airfoil=naca0008,
            ),
        ],
    )
    ### Build the fuselage geometry
    blend = lambda x: (1 - np.cos(np.pi * x)) / 2
    fuse_x_c = []
    fuse_z_c = []
    fuse_radius = []
    fuse_resolution = 10
    # Nose geometry
    fuse_nose_theta = np.linspace(0, np.pi / 2, fuse_resolution)
    fuse_x_c.extend(
        [
            (wing_x_quarter_chord - wing_root_chord / 4) - nose_length * np.cos(theta)
            for theta in fuse_nose_theta
        ]
    )
    fuse_z_c.extend([-fuse_diameter / 2] * fuse_resolution)
    fuse_radius.extend([fuse_diameter / 2 * np.sin(theta) for theta in fuse_nose_theta])
    # Taper
    fuse_taper_x_nondim = np.linspace(0, 1, fuse_resolution)
    fuse_x_c.extend(
        [
            0.0 * boom_length + (0.6 - 0.0) * boom_length * x_nd
            for x_nd in fuse_taper_x_nondim
        ]
    )
    fuse_z_c.extend(
        [
            -fuse_diameter / 2 * blend(1 - x_nd) - boom_diameter / 2 * blend(x_nd)
            for x_nd in fuse_taper_x_nondim
        ]
    )
    fuse_radius.extend(
        [
            fuse_diameter / 2 * blend(1 - x_nd) + boom_diameter / 2 * blend(x_nd)
            for x_nd in fuse_taper_x_nondim
        ]
    )
    # Tail
    # fuse_tail_x_nondim = np.linspace(0, 1, fuse_resolution)[1:]
    # fuse_x_c.extend([
    #     0.9 * boom_length + (1 - 0.9) * boom_length * x_nd for x_nd in fuse_taper_x_nondim
    # ])
    # fuse_z_c.extend([
    #     -boom_diameter / 2 * blend(1 - x_nd) for x_nd in fuse_taper_x_nondim
    # ])
    # fuse_radius.extend([
    #     boom_diameter / 2 * blend(1 - x_nd) for x_nd in fuse_taper_x_nondim
    # ])
    fuse_straight_resolution = 4
    fuse_x_c.extend(
        [
            0.6 * boom_length + (1 - 0.6) * boom_length * x_nd
            for x_nd in np.linspace(0, 1, fuse_straight_resolution)[1:]
        ]
    )
    fuse_z_c.extend([-boom_diameter / 2] * (fuse_straight_resolution - 1))
    fuse_radius.extend([boom_diameter / 2] * (fuse_straight_resolution - 1))

    fuse = asb.Fuselage(
        name="Fuselage",
        x_le=0,
        y_le=0,
        z_le=0,
        xsecs=[
            asb.FuselageXSec(x_c=fuse_x_c[i], z_c=fuse_z_c[i], radius=fuse_radius[i])
            for i in range(len(fuse_x_c))
        ],
    )

    # Assemble the airplane
    fuses = []
    hstabs = []
    vstabs = []
    if n_booms == 1:
        fuses.append(fuse)
        hstabs.append(hstab)
        vstabs.append(vstab)
    elif n_booms == 2:
        boom_location = 0.40  # as a fraction of the half-span

        left_fuse = copy.deepcopy(fuse)
        right_fuse = copy.deepcopy(fuse)
        left_fuse.xyz_le += cas.vertcat(0, -wing_span / 2 * boom_location, 0)
        right_fuse.xyz_le += cas.vertcat(0, wing_span / 2 * boom_location, 0)
        fuses.extend([left_fuse, right_fuse])

        left_hstab = copy.deepcopy(hstab)
        right_hstab = copy.deepcopy(hstab)
        left_hstab.xyz_le += cas.vertcat(0, -wing_span / 2 * boom_location, 0)
        right_hstab.xyz_le += cas.vertcat(0, wing_span / 2 * boom_location, 0)
        hstabs.extend([left_hstab, right_hstab])

        left_vstab = copy.deepcopy(vstab)
        right_vstab = copy.deepcopy(vstab)
        left_vstab.xyz_le += cas.vertcat(0, -wing_span / 2 * boom_location, 0)
        right_vstab.xyz_le += cas.vertcat(0, wing_span / 2 * boom_location, 0)
        vstabs.extend([left_vstab, right_vstab])

    elif n_booms == 3:
        boom_location = 0.57  # as a fraction of the half-span

        left_fuse = copy.deepcopy(fuse)
        center_fuse = copy.deepcopy(fuse)
        right_fuse = copy.deepcopy(fuse)
        left_fuse.xyz_le += cas.vertcat(0, -wing_span / 2 * boom_location, 0)
        right_fuse.xyz_le += cas.vertcat(0, wing_span / 2 * boom_location, 0)
        fuses.extend([left_fuse, center_fuse, right_fuse])

        left_hstab = copy.deepcopy(hstab)
        center_hstab = copy.deepcopy(hstab)
        right_hstab = copy.deepcopy(hstab)
        left_hstab.xyz_le += cas.vertcat(0, -wing_span / 2 * boom_location, 0)
        right_hstab.xyz_le += cas.vertcat(0, wing_span / 2 * boom_location, 0)
        hstabs.extend([left_hstab, center_hstab, right_hstab])

        left_vstab = copy.deepcopy(vstab)
        center_vstab = copy.deepcopy(vstab)
        right_vstab = copy.deepcopy(vstab)
        left_vstab.xyz_le += cas.vertcat(0, -wing_span / 2 * boom_location, 0)
        right_vstab.xyz_le += cas.vertcat(0, wing_span / 2 * boom_location, 0)
        vstabs.extend([left_vstab, center_vstab, right_vstab])

    else:
        raise ValueError("Bad value of n_booms!")

    airplane = asb.Airplane(
        name="Solar1",
        x_ref=0,
        y_ref=0,
        z_ref=0,
        wings=[wing] + hstabs + vstabs,
        fuselages=fuses,
    )

    return airplane
