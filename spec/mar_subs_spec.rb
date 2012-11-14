require_relative '../lib/mar_subs'

describe 'normalize_row' do
  specify do
  control = 2.0
  data = [1, 1.0, '1.0', -1, -1.0, '-1.0', nil]
  expect = [0.5, 0.5, 0.5, -0.5, -0.5, -0.5, nil]
  normalize_row(2.0, data).should eql expect
  normalize_row('2.0', data).should eql expect
  end
end
  