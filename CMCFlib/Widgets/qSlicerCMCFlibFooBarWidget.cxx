/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerCMCFlibFooBarWidget.h"
#include "ui_qSlicerCMCFlibFooBarWidget.h"

//-----------------------------------------------------------------------------
class qSlicerCMCFlibFooBarWidgetPrivate
  : public Ui_qSlicerCMCFlibFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerCMCFlibFooBarWidget);
protected:
  qSlicerCMCFlibFooBarWidget* const q_ptr;

public:
  qSlicerCMCFlibFooBarWidgetPrivate(
    qSlicerCMCFlibFooBarWidget& object);
  virtual void setupUi(qSlicerCMCFlibFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerCMCFlibFooBarWidgetPrivate
::qSlicerCMCFlibFooBarWidgetPrivate(
  qSlicerCMCFlibFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerCMCFlibFooBarWidgetPrivate
::setupUi(qSlicerCMCFlibFooBarWidget* widget)
{
  this->Ui_qSlicerCMCFlibFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerCMCFlibFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerCMCFlibFooBarWidget
::qSlicerCMCFlibFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerCMCFlibFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerCMCFlibFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerCMCFlibFooBarWidget
::~qSlicerCMCFlibFooBarWidget()
{
}
