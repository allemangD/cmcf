from __future__ import annotations

import time
from typing import TYPE_CHECKING, Any, Optional

import slicer
from MRMLCorePython import (
    vtkMRMLModelNode,
    vtkMRMLSequenceNode,
    vtkMRMLSubjectHierarchyNode,
    vtkMRMLScene, vtkMRMLModelDisplayNode,
)
from slicer.ScriptedLoadableModule import *
from slicer.i18n import tr as _
from slicer.i18n import translate
from slicer.util import VTKObservationMixin
from vtkSlicerSequencesModuleLogicPython import vtkSlicerSequencesLogic
from vtkSlicerSequencesModuleMRMLPython import vtkMRMLSequenceBrowserNode

if TYPE_CHECKING:
    from vtkSlicerCMCFlibModuleLogicPython import vtkSlicerCMCFLibLogic


class CMCF(ScriptedLoadableModule):
    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = _("CMCF")
        self.parent.categories = [translate("qSlicerAbstractCoreModule", "CMCF")]
        self.parent.dependencies = ["CMCFlib", "Sequences"]
        self.parent.contributors = [
            "David Allemang (University of North Carolina at Chapel Hill)"
        ]
        self.parent.helpText = _("")
        self.parent.acknowledgementText = _("")


class CMCFUI:
    # Qt type annotations don't work, so just note these with comments.

    selModel: Any  # qMRMLSubjectHierarchyTreeView
    dsbRate: Any  # QDoubleSpinBox
    sbStages: Any  # QSpinBox
    selSequence: Any  # qMRMLNodeComboBox
    btnApply: Any  # QPushButton
    btnIdentify: Any  # QPushButton


class CMCFWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    ui: CMCFUI
    logic: CMCFLogic
    lib_logic: vtkSlicerCMCFLibLogic

    _model: Optional[vtkMRMLModelNode] = None
    _sequence: Optional[vtkMRMLSequenceNode] = None

    def __init__(self, parent=None) -> None:
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)

    def setup(self) -> None:
        """Called when the user opens the module the first time and the widget is initialized."""
        ScriptedLoadableModuleWidget.setup(self)

        root_widget = slicer.util.loadUI(self.resourcePath("UI/CMCF.ui"))
        self.layout.addWidget(root_widget)
        self.ui = CMCFUI()
        self.ui.__dict__.update(
            {  # Hack so that completions work in Python interactor.
                widget.name: widget
                for widget in slicer.util.findChildren(root_widget)
                if getattr(widget, "name", None)
            }
        )

        root_widget.setMRMLScene(slicer.mrmlScene)

        self.logic = CMCFLogic()
        self.lib_logic = slicer.modules.cmcflib.logic()  # noqa

        # Connections
        self.ui.selModel.currentItemChanged.connect(self._onChangedModel)
        self.ui.selSequence.currentNodeChanged.connect(self._onChangedSequence)
        self.ui.btnApply.clicked.connect(self._onApply)
        self.ui.btnIdentify.clicked.connect(self._onIdentify)

    def cleanup(self) -> None:
        """Called when the application closes and the module widget is destroyed."""
        pass

    def enter(self) -> None:
        """Called each time the user opens this module."""
        pass

    def exit(self) -> None:
        """Called each time the user opens a different module."""
        pass

    def _onChangedModel(self, item: int):
        sh: vtkMRMLSubjectHierarchyNode = slicer.mrmlScene.GetSubjectHierarchyNode()

        if not item:
            self.ui.btnApply.enabled = False
            self.ui.selSequence.setCurrentNode(None)
            self._model = None
        else:
            node = sh.GetItemDataNode(item)
            self.ui.btnApply.enabled = True
            self.ui.selSequence.setCurrentNode(node.GetNodeReference("CMCF_SEQUENCE"))
            self._model = node

    def _onChangedSequence(self, node: vtkMRMLSequenceNode):
        self._sequence = node

        self.ui.btnIdentify.enabled = bool(self._sequence)

    def _onApply(self):
        """Invoke vtkSlicerCMCFLibLogic::GenerateCMCFSequence and update displays with the result."""

        if self._sequence is None:
            node = self._sequence = vtkMRMLSequenceNode()
            node.SetName(self._model.GetName() + " CMCF")
            slicer.mrmlScene.AddNode(node)
            self.ui.selSequence.setCurrentNode(node)

        self._model.SetNodeReferenceID(
            "CMCF_SEQUENCE", self._sequence and self._sequence.GetID()
        )

        with slicer.util.WaitCursor():
            _start = time.time()
            rate = self.ui.dsbRate.value
            stages = self.ui.sbStages.value
            self.lib_logic.GenerateCMCFSequence(
                self._model, self._sequence, rate, stages
            )
            _end = time.time()
            delta = int((_end - _start) * 1000)
            slicer.util.showStatusMessage(f"CMCF: Completed {stages} in {delta}ms.")

        browser: vtkMRMLSequenceBrowserNode = (
            slicer.modules.sequences.toolBar().activeBrowserNode()
        )
        if not browser:
            browser = slicer.mrmlScene.AddNode(vtkMRMLSequenceBrowserNode())
            slicer.modules.sequences.toolBar().setActiveBrowserNode(browser)

        try:
            slicer.mrmlScene.StartState(vtkMRMLScene.BatchProcessState)
            sequences: vtkSlicerSequencesLogic = slicer.modules.sequences.logic()
            sequences.AddSynchronizedNode(self._sequence, None, browser)

            proxy: vtkMRMLModelNode = browser.GetProxyNode(self._sequence)
            proxy.SaveWithSceneOff()
        finally:
            slicer.mrmlScene.EndState(vtkMRMLScene.BatchProcessState)

    def _onIdentify(self):
        if self._sequence is None:
            slicer.util.errorDisplay(
                "Missing sequence data. Select a model and compute a sequence first.",
            )
            return

        import builtins
        from vtkmodules.numpy_interface.dataset_adapter import WrapDataObject

        builtins._result = WrapDataObject(
            result := self.lib_logic.IdentifyParabolics(self._sequence),
        )

        slicer.mrmlScene.AddNode(model := vtkMRMLModelNode())
        model.SetAndObservePolyData(result)
        model.CreateDefaultDisplayNodes()
        disp: vtkMRMLModelDisplayNode = model.GetDisplayNode()
        disp.SetColor(0, 0, 0)
        disp.SetOpacity(0.8)
        # disp.SetRepresentation(vtkMRMLModelDisplayNode.PointsRepresentation)
        disp.SetRepresentation(vtkMRMLModelDisplayNode.WireframeRepresentation)


class CMCFLogic(ScriptedLoadableModuleLogic, VTKObservationMixin):
    def __init__(self, parent=None) -> None:
        ScriptedLoadableModuleLogic.__init__(self, parent)
        VTKObservationMixin.__init__(self)


class CMCFTest(ScriptedLoadableModuleTest):
    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
