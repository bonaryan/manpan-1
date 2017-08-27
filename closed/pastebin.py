##### Scripts to copy-paste in the parametric job run... #####

### get the radius for a lab specimen and produce a model in ABQ
radius = closed_polygons.class_2_radius(
    n_sides = parameter[0],
    p_classification = parameter[1],
    thickness = 3,
    f_yield = 420.
    )

job_return = closed_polygons.modeler(
    n_sides = parameter[0], 
    r_circle = radius, 
    p_classification = parameter[1], 
    n_of_waves = parameter[0] / 2, 
    m_of_waves = parameter[0] / 2, 
    u_max = 0.00667, 
    f_yield = 420., 
    IDstring = job_ID, 
    proj_name = prj_name, 
    )

### Calculate data for lab specimens (given thickness-calss, calculated radius)
import abq_toolset
import EN_tools
thickness = 3
r_circle = closed_polygons.class_2_radius(
    n_sides = parameter[0],
    p_classification = parameter[1],
    thickness = thickness,
    f_yield = 420.
    )

profile = closed_polygons.cs_calculator(
    n_sides = parameter[0], 
    r_circle = r_circle, 
    p_classification = parameter[1], 
    column_length = 1000.,
    f_yield = 420., 
    )

x_corners = profile[-2]
y_corners = profile[-1]

nodes = [x_corners, y_corners]
elem = [range(0, len(x_corners)), range(1, len(x_corners)) + [0], len(x_corners) * [thickness]]
cs_properties = abq_toolset.cs_prop(nodes, elem)
area = cs_properties[0]
I_2 = cs_properties[-2]

EN_tools.lmda = (
    1000.,
    area,
    I_2,
    kapa_BC = 1.,
    E_modulus = 210000.,
    f_yield = 420.
    )

job_return = [r_circle]+[profile[0:3]+[lmda]]



###### create and run models

job_return = closed_polygons.modeler(
        n_sides = parameter[0], 
        p_classification = parameter[1], 
        n_of_waves = parameter[0] / 2, 
        m_of_waves = parameter[0] / 2, 
        u_max = parameter[3] * 1000, 
        f_yield = 760, 
        nominal_fy = '650'
        IDstring = job_ID, 
        proj_name = prj_name, 
        )