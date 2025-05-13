import slicer
from slicer.ScriptedLoadableModule import *
from slicer.i18n import tr as _
from slicer.i18n import translate
from slicer.util import VTKObservationMixin


class CMCF(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("CMCF")
        self.parent.categories = [translate("qSlicerAbstractCoreModule", "CMCF")]
        self.parent.dependencies = []
        self.parent.contributors = ["David Allemang (UNC)"]
        self.parent.helpText = _("")
        self.parent.acknowledgementText = _("")


class CMCFWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    def __init__(self, parent=None) -> None:
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)

        # self.logic =
        self.logic = None

    def setup(self) -> None:
        """Called when the user opens the module the first time and the widget is initialized."""
        ScriptedLoadableModuleWidget.setup(self)

        # Load widget from .ui file (created by Qt Designer).
        # Additional widgets can be instantiated manually and added to self.layout.
        uiWidget = slicer.util.loadUI(self.resourcePath("UI/CMCF.ui"))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)

        # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
        # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
        # "setMRMLScene(vtkMRMLScene*)" slot.
        uiWidget.setMRMLScene(slicer.mrmlScene)

        # Create logic class. Logic implements all computations that should be possible to run
        # in batch mode, without a graphical user interface.
        self.logic = CMCFLogic()

        # Connections

        # # These connections ensure that we update parameter node when scene is closed
        # self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
        # self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

        # # Buttons
        # self.ui.applyButton.connect("clicked(bool)", self.onApplyButton)

    def cleanup(self) -> None:
        """Called when the application closes and the module widget is destroyed."""
        pass

    def enter(self) -> None:
        """Called each time the user opens this module."""
        pass

    def exit(self) -> None:
        """Called each time the user opens a different module."""
        pass

    def onSceneStartClose(self, caller, event) -> None:
        """Called just before the scene is closed."""
        pass

    def onSceneEndClose(self, caller, event) -> None:
        """Called just after the scene is closed."""
        pass


class CMCFLogic(ScriptedLoadableModuleLogic, VTKObservationMixin):
    def __init__(self, parent=None) -> None:
        ScriptedLoadableModuleLogic.__init__(self, parent)
        VTKObservationMixin.__init__(self)

    def getParameterNode(self):
        # return CMCFParameterNode(super().getParameterNode())
        return super().getParameterNode()


class CMCFTest(ScriptedLoadableModuleTest):
    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
