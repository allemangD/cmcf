from MRMLCorePython import vtkMRMLModelNode, vtkMRMLSequenceNode


class vtkSlicerCMCFLibLogic:
    @staticmethod
    def GenerateCMCFSequence(
            model: vtkMRMLModelNode,
            sequence: vtkMRMLSequenceNode,
            rate: float = 0.05,
            stages: int = 100,
    ):
        ...

    @staticmethod
    def IdentifyParabolics(
            sequence: vtkMRMLSequenceNode,
            skip: int = 5,
            tolerance: float = 0.005,
    ):
        ...
