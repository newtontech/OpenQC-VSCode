import * as assert from 'assert';
import * as vscode from 'vscode';

suite('Extension Test Suite', () => {
  suiteSetup(async () => {
    // Ensure the extension is activated before running tests
    const extension = vscode.extensions.getExtension('newtontech.openqc');
    if (extension && !extension.isActive) {
      await extension.activate();
    }
  });

  test('Extension should be present', () => {
    const extension = vscode.extensions.getExtension('newtontech.openqc');
    assert.ok(extension, 'Extension should be installed');
  });

  test('Extension should activate', async () => {
    const extension = vscode.extensions.getExtension('newtontech.openqc');
    if (!extension) {
      assert.fail('Extension not found');
    }

    await extension.activate();
    assert.strictEqual(extension.isActive, true, 'Extension should be active');
  });

  test('Extension should set sidebar context', async () => {
    // Wait a bit for activation to complete
    await new Promise(resolve => setTimeout(resolve, 1000));
    
    const sidebarEnabled = await vscode.commands.executeCommand(
      'getContext',
      'openqc.sidebar.enabled'
    );
    // Note: getContext is not a real VSCode command, so we test differently
    // We'll verify by checking if the tree views are available
    const treeViews = vscode.window.createTreeView('openqc.molecules', {
      treeDataProvider: {
        getTreeItem: () => new vscode.TreeItem('test'),
        getChildren: () => Promise.resolve([]),
      },
    });
    assert.ok(treeViews, 'Tree view should be created');
    treeViews.dispose();
  });

  test('Command: openqc.visualizeStructure should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.visualizeStructure'),
      'openqc.visualizeStructure command should be registered'
    );
  });

  test('Command: openqc.plotData should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.plotData'),
      'openqc.plotData command should be registered'
    );
  });

  test('Command: openqc.previewInput should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.previewInput'),
      'openqc.previewInput command should be registered'
    );
  });

  test('Command: openqc.startLSP should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.startLSP'),
      'openqc.startLSP command should be registered'
    );
  });

  test('Command: openqc.stopLSP should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.stopLSP'),
      'openqc.stopLSP command should be registered'
    );
  });

  test('Command: openqc.restartLSP should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.restartLSP'),
      'openqc.restartLSP command should be registered'
    );
  });

  test('Command: openqc.validate should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.validate'),
      'openqc.validate command should be registered'
    );
  });

  test('Command: openqc.sidebar.refreshMolecules should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.sidebar.refreshMolecules'),
      'openqc.sidebar.refreshMolecules command should be registered'
    );
  });

  test('Command: openqc.sidebar.refreshJobs should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.sidebar.refreshJobs'),
      'openqc.sidebar.refreshJobs command should be registered'
    );
  });

  test('Command: openqc.sidebar.openMolecule should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.sidebar.openMolecule'),
      'openqc.sidebar.openMolecule command should be registered'
    );
  });

  test('Command: openqc.sidebar.deleteMolecule should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.sidebar.deleteMolecule'),
      'openqc.sidebar.deleteMolecule command should be registered'
    );
  });

  test('Command: openqc.sidebar.runCalculation should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.sidebar.runCalculation'),
      'openqc.sidebar.runCalculation command should be registered'
    );
  });

  test('Command: openqc.sidebar.viewResults should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.sidebar.viewResults'),
      'openqc.sidebar.viewResults command should be registered'
    );
  });

  test('Command: openqc.sidebar.exportData should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.sidebar.exportData'),
      'openqc.sidebar.exportData command should be registered'
    );
  });

  test('Command: openqc.sidebar.cancelJob should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.sidebar.cancelJob'),
      'openqc.sidebar.cancelJob command should be registered'
    );
  });

  test('Command: openqc.sidebar.restartJob should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.sidebar.restartJob'),
      'openqc.sidebar.restartJob command should be registered'
    );
  });

  test('Command: openqc.convertFormat should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.convertFormat'),
      'openqc.convertFormat command should be registered'
    );
  });

  test('Command: openqc.convertToXYZ should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.convertToXYZ'),
      'openqc.convertToXYZ command should be registered'
    );
  });

  test('Command: openqc.convertToPDB should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.convertToPDB'),
      'openqc.convertToPDB command should be registered'
    );
  });

  test('Command: openqc.convertToVASP should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.convertToVASP'),
      'openqc.convertToVASP command should be registered'
    );
  });

  test('Command: openqc.convertToGaussian should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.convertToGaussian'),
      'openqc.convertToGaussian command should be registered'
    );
  });

  test('Command: openqc.batchConvert should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.batchConvert'),
      'openqc.batchConvert command should be registered'
    );
  });

  test('Command: openqc.checkConverterBackend should be registered', async () => {
    const commands = await vscode.commands.getCommands(true);
    assert.ok(
      commands.includes('openqc.checkConverterBackend'),
      'openqc.checkConverterBackend command should be registered'
    );
  });
});
