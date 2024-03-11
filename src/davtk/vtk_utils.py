import vtk

def new_prop(ns):
    prop = vtk.vtkProperty()
    update_prop(prop, ns)
    return prop

def update_prop(prop, ns):
    try:
        if ns.color is not None:
            prop.SetColor(ns.color)
    except AttributeError:
        pass

    try:
        if ns.opacity is not None:
            prop.SetOpacity(ns.opacity)
    except AttributeError:
        pass

    try:
        if ns.specular is not None:
            prop.SetSpecular(ns.specular)
    except AttributeError:
        pass

    try:
        if ns.specular_radius is not None:
            prop.SetSpecularPower(1.0/ns.specular_radius)
    except AttributeError:
        pass

    try:
        if ns.ambient is not None:
            prop.SetAmbient(ns.ambient)
    except AttributeError:
        pass
