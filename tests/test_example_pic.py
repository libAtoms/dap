import os
import shutil
import subprocess
import vtk
from pathlib import Path

def test_example(tmp_path):
    os.chdir(tmp_path)
    assets = Path(__file__).parent / 'assets'
    for f in ["example.xyz", "example.settings", "example.commands", "example.iso"]:
        shutil.copy(assets / f, f)

    ref_reader = vtk.vtkPNGReader()
    ref_reader.SetFileName(str(assets / "example.png"))
    ref_reader.Update()
    ref_img = vtk.vtkImageExtractComponents()
    ref_img.SetComponents(0,1,2)
    ref_img.SetInputConnection(ref_reader.GetOutputPort())
    ref_img.Update()

    p = subprocess.Popen(['dap', '-e', 'read example.settings', 
                             '-e', 'read example.commands',
                             '-e', 'view -dir 1 0.1 0.2  0 0 1',
                             'example.xyz'],
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate('snapshot test.png\nexit\n'.encode())
    assert p.returncode == 0

    test_reader = vtk.vtkPNGReader()
    test_reader.SetFileName("test.png")
    test_reader.Update()
    test_img = vtk.vtkImageExtractComponents()
    test_img.SetComponents(0,1,2)
    test_img.SetInputConnection(test_reader.GetOutputPort())
    test_img.Update()

    idiff = vtk.vtkImageDifference()
    # idiff.SetThreshold(1)
    idiff.AllowShiftOff()
    idiff.SetImageData(ref_img.GetOutput())
    idiff.SetInputConnection(test_img.GetOutputPort())
    idiff.Update()

    # print("thresh error", idiff.GetThresholdedError())
    assert idiff.GetError() < 1.0e-10
