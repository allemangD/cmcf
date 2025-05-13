/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// Qt includes
#include <QDebug>

// Slicer includes
#include "qSlicerCMCFlibModuleWidget.h"
#include "ui_qSlicerCMCFlibModuleWidget.h"

//-----------------------------------------------------------------------------
class qSlicerCMCFlibModuleWidgetPrivate: public Ui_qSlicerCMCFlibModuleWidget
{
public:
  qSlicerCMCFlibModuleWidgetPrivate();
};

//-----------------------------------------------------------------------------
// qSlicerCMCFlibModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerCMCFlibModuleWidgetPrivate::qSlicerCMCFlibModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerCMCFlibModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerCMCFlibModuleWidget::qSlicerCMCFlibModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerCMCFlibModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerCMCFlibModuleWidget::~qSlicerCMCFlibModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerCMCFlibModuleWidget::setup()
{
  Q_D(qSlicerCMCFlibModuleWidget);
  d->setupUi(this);
  this->Superclass::setup();
}
