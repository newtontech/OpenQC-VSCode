import { renderResultsWebviewHtml } from '../../../src/webviews/resultsWebview';

describe('resultsWebview', () => {
  it('renders the results template with a strict CSP and nonce-protected styles', () => {
    const html = renderResultsWebviewHtml(
      {
        jobId: 'job-1',
        jobName: 'Water Optimization',
        software: 'Gaussian',
        status: 'completed',
        startTime: '2026-06-07T10:00:00.000Z',
        endTime: '2026-06-07T10:05:00.000Z',
        duration: '5 minutes',
        output: 'Final energy: -76.0',
      },
      'vscode-resource:'
    );

    const nonceMatch = html.match(/<style nonce="([A-Za-z0-9]{32})">/);

    expect(html).toContain('<meta http-equiv="Content-Security-Policy"');
    expect(html).toContain("default-src 'none';");
    expect(html).toContain("script-src 'none';");
    expect(html).toContain('img-src vscode-resource: data:;');
    expect(html).toContain('font-src vscode-resource:;');
    expect(nonceMatch).not.toBeNull();
    expect(html).toContain(`style-src 'nonce-${nonceMatch?.[1]}';`);
    expect(html).toContain('background: #3c9;');
    expect(html).toContain('<span class="info-label">Started:</span>');
    expect(html).toContain('<span class="info-label">Completed:</span>');
  });

  it('escapes rendered job data before interpolating it into the template', () => {
    const html = renderResultsWebviewHtml(
      {
        jobId: `job-"'><img src=x onerror=alert(1)>`,
        jobName: '<script>alert("name")</script>',
        software: 'Gaussian & ORCA',
        status: 'failed',
        duration: `<b>instant</b>`,
        output: `line 1\n<script>alert('output')</script>`,
      },
      'vscode-resource:'
    );

    expect(html).toContain('&lt;script&gt;alert(&quot;name&quot;)&lt;/script&gt;');
    expect(html).toContain('job-&quot;&#39;&gt;&lt;img src=x onerror=alert(1)&gt;');
    expect(html).toContain('Gaussian &amp; ORCA');
    expect(html).toContain('&lt;b&gt;instant&lt;/b&gt;');
    expect(html).toContain('&lt;script&gt;alert(&#39;output&#39;)&lt;/script&gt;');
    expect(html).not.toContain('<script>alert');
    expect(html).not.toContain('<img src=x');
    expect(html).toContain('background: #f66;');
  });

  it('omits optional timestamps and defaults missing duration', () => {
    const html = renderResultsWebviewHtml(
      {
        jobId: 'job-2',
        jobName: 'Queued Job',
        software: 'CP2K',
        status: 'failed',
        output: 'No duration was captured',
      },
      'vscode-resource:'
    );

    expect(html).toContain('<span class="info-value">N/A</span>');
    expect(html).not.toContain('<span class="info-label">Started:</span>');
    expect(html).not.toContain('<span class="info-label">Completed:</span>');
  });
});
