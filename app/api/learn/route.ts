import { NextResponse } from 'next/server'
import { getLearnResources } from '@/lib/server-markdown-loader'

export async function GET() {
  try {
    const resources = getLearnResources()
    return NextResponse.json(resources)
  } catch (error) {
    console.error('Error fetching learn resources:', error)
    return NextResponse.json({ error: 'Failed to fetch learn resources' }, { status: 500 })
  }
}